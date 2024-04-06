import csv
from pathlib import Path
from typing import Dict, Any, Callable, List
from io import StringIO
from frozendict import frozendict
import warnings
import pint
from numpy import exp, power
from itertools import combinations
from uncertainties import ufloat, unumpy
import matplotlib.pyplot as plt
import numpy as np

""" This is a module for reading and calculating thermophysical properties from the MSTDB using ideal estimations (e.g. 
additivity of molar volumes) and the RK expansion for estimating the effects nonideal mixing.

Author: Matthew Louis
Email: matthewlouis31@gmail.com
"""

ureg = pint.UnitRegistry(auto_reduce_dimensions=True)
ureg.default_format = "~P"
Q_ = ureg.Quantity

# ---------------------
# Customized Warnings
# ---------------------

# ANSI escape code for yellow
YELLOW = '\033[93m'
RESET = '\033[0m'

# Customize warning format and add color
warnings.formatwarning = lambda message, category, filename, lineno, line=None: \
    f'{YELLOW}{filename}:{lineno}: {category.__name__}: {message}{RESET}\n'


def docstring_injector(cls):
    """This dectorator iterates overall all the methods of a class looking for callable objects
    with docstrings. It then replaces the docstring with a formatted string that includes the unique tp names."""
    unique_tp_names = ', '.join(cls._UNIQUE_TP_NAMES)
    for name, method in cls.__dict__.items():
        if callable(method) and method.__doc__:
            method.__doc__ = method.__doc__.format(UNIQUE_TP_NAMES=unique_tp_names)
    return cls

class ArbitraryThermoFunction:
    """This class is used for defining functions that represent the thermophysical properties of molten salts. The objects of this
    class are functions of a float (temperature) and they output the desired thermophysical property. This class conveniently
    defines addition, multiplication, etc. and performs error propagation to ensure that the result ThermoFunction has the correct
    applicable temperature range (or uncertainty)."""

    def __init__(self, func, tmin=None, tmax=None, uncertainty=None, units=None):
        self.func = func
        self.min_temp = tmin
        self.max_temp = tmax
        self.fractional_uncertainty = uncertainty
        self.units = units

    def __call__(self, temp):
        if self.min_temp is not None and self.max_temp is not None:
            if temp < self.min_temp or temp > self.max_temp:
                warnings.warn(
                    f"Temperature {temp} is outside the valid range [{self.min_temp}, {self.max_temp}]"
                , UserWarning)
        result = self.func(temp)
        unc = abs(result*self.fractional_uncertainty) if self.fractional_uncertainty is not None else 0.0
        try:
            return ufloat(result, unc)
        except: # Result is already a ufloat or another non-float datatype
            try: # When using np.sum(), results can get converted into a np array of ufloats, just use .item() to conver to ufloat
                return result.item()
            except: # Result is already a ufloat
                return result

    def __add__(self, other):
        if isinstance(other, ArbitraryThermoFunction):
            new_min_temp = max(self.min_temp, other.min_temp) if self.min_temp is not None and other.min_temp is not None else None
            new_max_temp = min(self.max_temp, other.max_temp) if self.max_temp is not None and other.max_temp is not None else None
            
            def new_func(temp):
                return self(temp) + other(temp)
            
            return ArbitraryThermoFunction(new_func, new_min_temp, new_max_temp)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            def new_func(temp):
                return self(temp) * other
            
            return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            def new_func(temp):
                return self(temp) / other
            
            return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            def new_func(temp):
                return other / self(temp)
            
            return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)
        else:
            return NotImplemented

    def reciprocal(self):
        def new_func(temp):
            return 1 / self(temp)
        
        return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)

    def log(self):
        def new_func(temp):
            return unumpy.log(self(temp))
        
        return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)

    def exp(self):
        def new_func(temp):
            return unumpy.exp(self(temp))
        
        return ArbitraryThermoFunction(new_func, self.min_temp, self.max_temp)

class ThermoFunction(ArbitraryThermoFunction):
    # Once the database is processed once, there is no longer a 'viscosity_exp' or 'viscosity_base10' key, instead they have just been
    # replaced with 'viscosity'. This class attribute keeps track of when this occurs.
    has_been_parsed = False
    
    def __init__(self, salt, key):
        self.salt = salt
        self.key = key
        
        # If dealing with viscosity, since the database is modified to exclude one of the viscosity keys
        # and have only 'viscosity', we need to catch this case
        if 'viscosity' in self.key and ThermoFunction.has_been_parsed:
            key_of_thermo_func = 'viscosity'
        else:
            key_of_thermo_func = self.key
        
        if isinstance(salt[key_of_thermo_func], ArbitraryThermoFunction):
            # If salt[key_of_thermo_func] is already a ThermoFunction, no need to parse
            # they min/max temp, etc. just take it from the ThermoFunction attributes
            func = salt[key_of_thermo_func].func
            min_temp = salt[key_of_thermo_func].min_temp
            max_temp = salt[key_of_thermo_func].max_temp
            fractional_uncertainty = salt[key_of_thermo_func].fractional_uncertainty
            units = salt[key_of_thermo_func].units
        else:
            coef_array = salt[key][0]
            
            if isinstance(salt[key][1], tuple):
                if len(salt[key]) == 3:
                    fractional_uncertainty = salt[key][2]/100
                else:
                    fractional_uncertainty = None
                min_temp, max_temp = salt[key][1]
            elif isinstance(salt[key][1], float):
                if len(salt[key]) == 3: # For some reason, only one point has been given for the temp range
                    fractional_uncertainty = salt[key][2]/100
                else:
                    fractional_uncertainty = salt[key][1]/100
                min_temp = None
                max_temp = None
            else:
                fractional_uncertainty = None
                min_temp = None
                max_temp = None
            
            units = Database._TP_UNITS[key]
            
            if key == 'liquid_heat_capacity':
                MW = salt['molecular_weight']
                func = lambda temp: Database._TP_FUNCTIONS[key](MW, *coef_array, temp)
            else:
                func = lambda temp: Database._TP_FUNCTIONS[key](*coef_array, temp)
        
        super().__init__(func, min_temp, max_temp, fractional_uncertainty, units)

@docstring_injector
class Database:
    # Configuration for CSV headers
    _CSV_HEADERS = {
        'formula': 'Formula',
        'molecular_weight': 'Molecular Weight',
        'composition': 'Composition (Mole %)',
        'uncertainty': 'Uncertainty (%)',
        'temp_uncertainty': 'Uncertainty (K)',
        'references': 'References',
        'temp_range': 'Applicable Temperature range (K)',
        'melting_temp': 'Melting T (K)',
        'boiling_temp': 'Boiling T (K)',
        'number': '#',
        'density': 'Density (g/cm3):   A - BT(K)',
        'viscosity_exp': 'Viscosity (mN*s/m2): A*exp(B/(R*T(K)))',
        'viscosity_base10': 'Viscosity (mN*s/m2): 10^(A + B/T + C/T**2)',
        'thermal_conductivity': 'Thermal Conductivity (W/m K):  A + B*T(K)',
        'liquid_heat_capacity': "Heat Capacity of Liquid (J/K mol):        \n A + B*T(K) + C*T-2(K) + D*T2(K)"
    }
    # Physical constants
    R = 8.314 # J/(K*mol)

    # ------------------------------------------------------------------------------------------
    # These are the functional expansions of the thermophysical properties given in the MSTDB,
    # they are all constant, and dependent only on the version/format of the MSTDB-TP
    # ------------------------------------------------------------------------------------------

    @staticmethod
    def VISCOSITY_EXP(A, B, T):
        # Note this returns dynamic viscosity in units of Pa*s
        return ( A*exp(B/(Database.R*T)) )*ureg('mN*s/m^2').to('Pa*s').magnitude
    
    @staticmethod
    def VISCOSITY_BASE_10(A, B, C, T):
        # Note this returns dynamic viscosity in units of Pa*s
        return ( power( 10, A + B/T + C/T**2 ) )*ureg('mN*s/m^2').to('Pa*s').magnitude
    
    @staticmethod
    def DENSITY(A, B, T):
        # Returns the density in SI units
        return ( A - B*T )*ureg('g/cm^3').to('kg/m^3').magnitude
    
    @staticmethod
    def THERMAL_CONDUCTIVITY(A, B, T):
        # Returns thermal conducitivty in SI units
        return A + B*T
    
    @staticmethod
    def HEAT_CAPACITY(MW, A, B, C, D, T):
        # Returns heat capacity in SI units, i.e. J/(kg * K), hence why we have to divide
        # by the molecular weight in g/mol
        return ( ( A + B*T + C/T**2 + D*T**2 )/(MW*ureg('g/mol').to('kg/mol')) ).magnitude
    
    _TP_FUNCTIONS = {
        'viscosity_exp': VISCOSITY_EXP,
        'viscosity_base10': VISCOSITY_BASE_10,
        'density': DENSITY,
        'thermal_conductivity': THERMAL_CONDUCTIVITY,
        'liquid_heat_capacity': HEAT_CAPACITY
    }

    # ------------------------------------------------------------------------------------------
    # These are the functional expansions of the thermophysical properties of higher order salt
    # systems using the ideal mixing assumptions
    # ------------------------------------------------------------------------------------------
    
    @staticmethod
    def IDEAL_VISCOSITY(tp_of_endmembers: List[ThermoFunction], composition_endmembers: list) -> ThermoFunction:
        """Calculate the ideal viscosity function.

        The ideal viscosity is calculated using the logarithmic mixing rule:

            \\log \\mu_id = \\sum_(i=1)^n \\log(\\mu_i)*x_i

        where \\mu_id is the ideal viscosity, x_i is the mole fraction of
        the ith component, and \\mu_i is the viscosity of the ith component."""

        log_viscosities = [x_i * mu_i.log() for mu_i, x_i in zip(tp_of_endmembers, composition_endmembers)]
        return np.sum(log_viscosities).exp()

    @staticmethod
    def IDEAL_DENSITY(tp_of_endmembers: List[ThermoFunction], composition_endmembers: list, mw_endmembers: list) -> ThermoFunction:
        """Calculate the ideal density function.

        The ideal density is calculated using the ideal mixing rule:

            \\rho_id = \\sum_(i=1)^n x_i*M_i / ( \\sum_(i=1)^n x_i*M_i/\\rho_i )

        where \\rho_id is the ideal density, x_i is the mole fraction of the ith component,
        M_i is the molecular weight of each component (in kg/mol), and \\rho_i is the density of the ith component."""

        weighted_densities = [x_i * (M_i / rho_i) for rho_i, x_i, M_i in zip(tp_of_endmembers, composition_endmembers, mw_endmembers)]
        summed_weighted_densities = np.sum(weighted_densities)
        total_mass_fraction = np.sum([x_i * M_i for x_i, M_i in zip(composition_endmembers, mw_endmembers)])
        return total_mass_fraction / summed_weighted_densities

    @staticmethod
    def IDEAL_THERMAL_CONDUCTIVITY(tp_of_endmembers: List[ThermoFunction], composition_endmembers: list) -> ThermoFunction:
        """Calculate the ideal thermal conductivity function.

        The ideal thermal conductivity is calculated using the following rule:

            k_id = \\sum_(i=1)^n x_i*k_i

        where k_id is the ideal thermal conductivity, x_i is the mole fraction of
        the ith component, and k_i is the thermal conductivity of the ith component."""
        return np.sum([x_i * k_i for k_i, x_i in zip(tp_of_endmembers, composition_endmembers)])

    @staticmethod
    def IDEAL_HEAT_CAPACITY(tp_of_endmembers: List[ThermoFunction], composition_endmembers: list) -> ThermoFunction:
        """Calculate the ideal heat capacity for a given temperature.

        The ideal heat capacity is calculated using the following rule:

            c_id = \\sum_(i=1)^n x_i*c_i

        where c_id is the ideal heat capacity, x_i is the mole fraction of
        the ith component, and c_i is the heat capacity of the ith component."""

        return np.sum([x_i * c_i for c_i, x_i in zip(tp_of_endmembers, composition_endmembers)])

    IDEAL_FUNCTIONS = {
        'viscosity': IDEAL_VISCOSITY,
        'density': IDEAL_DENSITY,
        'thermal_conductivity': IDEAL_THERMAL_CONDUCTIVITY,
        'liquid_heat_capacity': IDEAL_HEAT_CAPACITY
    }

    _TP_UNITS = {
        'viscosity_exp': 'Pa*s',
        'viscosity_base10': 'Pa*s',
        'density': 'kg/m^3',
        'thermal_conductivity': 'W/(m*K)',
        'liquid_heat_capacity': 'J/(kg*K)'
    }
    
    _UNIQUE_TP_NAMES = [
        'viscosity',
        'density',
        'thermal_conductivity',
        'liquid_heat_capacity'
    ]

    _CSV_VALUES = {
        'pure_salt': 'Pure Salt'
    }
    # Headers for the RK expansion coefficient file
    _CSV_HEADERS_RK = {
        'component1': 'C1',
        'component2': 'C2',
        'a1': 'A1',
        'b1': 'B1',
        'a2': 'A2',
        'b2': 'B2',
        'a3': 'A3',
        'b3': 'B3',
        'tmin': 'T min',
        'tmax': 'T max',
        'reference': 'Reference'
    }

    # Composition tolerance for selecting higher order system directly from the database (without performing property estimation)
    COMPOSITION_TOLERANCE = 0.05

    class _Parsers:
        """This class contains all of the parsers that are necessary for processing the .csv files. It's easier
        to logically group these under a helper class than to clutter the body of the class."""

        def __init__(self, outer_instance):
            self._outer = outer_instance
            pass

        def _parse_mstdb_tp(self, mstdb_tp_path: Path) -> Dict:
            """Parse the MSTDB TP file and return the data as a dictionary with frozendict keys."""
            data = dict()
            with mstdb_tp_path.open('r') as csvfile:
                # Now read the rest of the csv
                reader = list(csv.DictReader(csvfile, delimiter=','))  # Assuming comma-separated values
                reader = self._preprocess_subheader(reader, csvfile)

                for row in reader:
                    composition_dict = self._parse_composition(row)
                    parsed_row = self._parse_row(row)
                    # Now, for each of the thermophysical properties with expansion functions, replace the array with a function
                    for tp_key in [ key for key in Database._TP_FUNCTIONS.keys() if key in parsed_row.keys() ]:
                        if parsed_row[tp_key] is not None: # Can be none if tp_key is the wrong viscosity functionalization
                            parsed_row[tp_key] = ThermoFunction(parsed_row, tp_key)
                        if 'viscosity' in tp_key:
                            parsed_row['viscosity'] = parsed_row.pop(tp_key)
                    data.update({composition_dict: parsed_row})
                ThermoFunction.has_been_parsed = True
            return data
        
        def _parse_mstdb_tp_rk(self, mstdb_tp_rk_path):
            """Function for parsing the .csv file containing the RK coefficients for the density expansion. Since this
            csv does not have a subheader, its processing should be much easier"""
            rk = dict()
            with mstdb_tp_rk_path.open('r', errors='replace') as csvfile:
                reader = csv.DictReader(csvfile, delimiter=',')
                for row in reader:
                    salt_pair_key = frozenset({row[Database._CSV_HEADERS_RK['component1']],\
                                            row[Database._CSV_HEADERS_RK['component2']]})
                    parsed_row = self._parse_row_rk(row)
                    rk.update( { salt_pair_key: parsed_row } )
                    
            return rk
        
        def _parse_row_rk(self, row):
            """Function for parsing rows of the MSTDB-TP RK csv"""
            parsed_row = [[],]
            keys_by_coefficient_order = [ ( Database._CSV_HEADERS_RK[f'a{i}'], Database._CSV_HEADERS_RK[f'b{i}'] ) for i in range(1,4) ]

            # First add coefficients for each order
            for keys_at_order in keys_by_coefficient_order:
                parsed_row[0].append([float(row[keys_at_order[0]]), float(row[keys_at_order[1]])])

            # Then add temperature range
            parsed_row.append((float(row[Database._CSV_HEADERS_RK['tmin']]), float(row[Database._CSV_HEADERS_RK['tmax']])))

            return parsed_row

        def _preprocess_subheader(self, reader, csvfile):
            """Function for correctly processing the subheaders in MSTDB-TP corresponding to the coefficients of the temperature
            dependence expansion of TP's"""
            # First remove the first row of reader, it corresponds to the subheader
            reader.pop(0)
            reader = [{key: value for key, value in row.items() if key != ''} for row in reader] # Filter out empty keys

            csvfile.seek(0) # Reset the cursor of the csv

            # Currently, since the Heat capacity of liquid header has an explicit newline (bad database design), it gets spilt
            # into two lines, so we have to handle this accordingly and append the headers
            headers = next(csvfile).strip().replace('\"', '').split(',')
            next_row = next(csvfile).strip().replace('\"', '').split(',')
            headers[-1] += '        \n ' + next_row[0]
            headers += next_row[1:]
            # Now, store the second row specially
            self._second_row = next(csvfile).strip().split(',')

            csvfile.seek(0) # Reset the cursor of the csv
            lines = csvfile.readlines()
            for row_index, row in enumerate(lines):
                if row_index in [0,1,2]:
                    continue # Skip header and subheader
                index_of_last_empty_col = -99
                for col_index, col in enumerate(self._smart_split(row.strip(), delimiter=',')):
                    if headers[col_index] == '': # We're in the subheader
                        index_of_last_empty_col = col_index
                        # So the value at the key of the last nonempty column should actually be a list
                        # Note we subtract by 3 because the first three rows correspond to the header and subheader
                        reader_at_header = reader[row_index-3][key_of_last_nonempty_col] 
                        if not isinstance(reader_at_header, list):
                            if reader_at_header == '----' or reader_at_header == '':
                                parsed_val = None
                            else:
                                parsed_val = float(reader_at_header)
                            reader[row_index-3][key_of_last_nonempty_col] = [[parsed_val], '']
                        # Now add this new value
                        if col == '----' or col == '':
                            append_val = None
                        else:
                            append_val = float(col)
                        reader[row_index-3][key_of_last_nonempty_col][0].append(append_val)
                    else:
                        if index_of_last_empty_col == col_index - 1:
                            # Add temp range and uncertainty to the end of the list
                            last_tp_is_viscosity =  key_of_last_nonempty_col in [Database._CSV_HEADERS['viscosity_exp'], \
                                                                                 Database._CSV_HEADERS['viscosity_base10']]
                            last_tp_coef_array = reader[row_index -3][key_of_last_nonempty_col][0]
                            if not last_tp_is_viscosity:
                                # If the list is all None's, just replace it with None
                                is_none = [coef == None for coef in last_tp_coef_array]
                                if all(is_none):
                                    reader[row_index -3][key_of_last_nonempty_col] = None
                                else:
                                    # Otherwise add on the temperature range and convert nones to zeros in the array
                                    temp_range = self._parse_temp_range(col)
                                    uncertainty = self._parse_temp_range(self._smart_split(lines[row_index].strip(), delimiter=',')[col_index + 1])
                                
                                    reader[row_index -3][key_of_last_nonempty_col][0] = [ 0 if coef is None else coef \
                                                                                         for coef in last_tp_coef_array ]
                                    reader[row_index -3][key_of_last_nonempty_col][1] = temp_range
                                    if uncertainty is not None:
                                        reader[row_index -3][key_of_last_nonempty_col].append(uncertainty)
                            else: # Last tp is viscosity
                                temp_range = self._parse_temp_range(col)
                                uncertainty = self._parse_temp_range(self._smart_split(lines[row_index].strip(), delimiter=',')[col_index + 1])

                                reader[row_index -3][key_of_last_nonempty_col][1] = temp_range
                                if uncertainty is not None:
                                    reader[row_index -3][key_of_last_nonempty_col].append(uncertainty)

                        key_of_last_nonempty_col = headers[col_index]

            csvfile.seek(0) # reset the cursor of the csv

            return reader

        def _smart_split(self, input_string, delimiter=',', quotechar='"'):
            """
            Splits a string by a delimiter, treating everything within quotechar as a single value.
            
            :param input_string: The string to split.
            :param delimiter: The delimiter to split by, defaults to ','.
            :param quotechar: The character to use as a quotation mark, defaults to '"'.
            :return: A list of values split by the delimiter, with quoted sections treated as single values.
            """
            # Use StringIO to turn the string into a file-like object
            stringio_obj = StringIO(input_string)
            
            # Create a CSV reader object
            reader = csv.reader(stringio_obj, delimiter=delimiter, quotechar=quotechar)
            
            # Since there's only one line in this string, we can use next() to get the data
            parsed_list = next(reader)
            
            return parsed_list

        def _parse_composition(self, row: dict):
            """Parse the composition string into a frozendict."""
            if row[Database._CSV_HEADERS['composition']] == Database._CSV_VALUES['pure_salt']:
                # The composition of the salt is just 100% of the listed formula for pure salts
                return frozendict({row[Database._CSV_HEADERS['formula']]: 1.0})
            else:
                composition = row[Database._CSV_HEADERS['composition']].split('-')
                formula = row[Database._CSV_HEADERS['formula']]
                endmembers = formula.split('-')
                composition_dict = {endmember: float(mole_frac_endmember) for endmember, mole_frac_endmember \
                                    in zip(endmembers, composition) }
                return frozendict(composition_dict)

        def _parse_row(self, row: Dict[str, str]) -> Dict[str, Any]:
            """Converts all numerical strings in the row to appropriate data types and handles missing values."""
            parsed_row = dict()
            for key, value in row.items():
                # Columns to skip, not currently relevant to output or parsed separately
                cols_to_skip = ['composition', 'temp_uncertainty', 'references', 'temp_range', 'number', 'formula']
                if key in [Database._CSV_HEADERS[col_header] for col_header in cols_to_skip  ]:  # Skip the composition column as it is handled separately
                    continue
                if value == '----':  # Placeholder for missing values
                    parsed_row[key] = None
                else:
                    try:
                        # Attempt to convert numerical values
                        parsed_row[key] = float(value)
                    except:
                        # Keep as string if conversion fails
                        parsed_row[key] = value
            # Now, since there are two possible expansions for viscosity, and data is given for only one, the other is None, and should be
            # excluded
            viscosity_keys = ['viscosity_exp', 'viscosity_base10']
            self._parse_viscosity_coefs(viscosity_keys, parsed_row)
            # Now convert all keys to the convenient names
            def get_first_key_with_value(value_to_match, my_dict):
                return next( ( key for key, value in my_dict.items() if value ==  value_to_match), None )
            parsed_row = {get_first_key_with_value(key, Database._CSV_HEADERS): value for key, value in parsed_row.items() }
            return parsed_row
        
        def _parse_viscosity_coefs(self, viscosity_keys, row):
            """This function is meant for processing viscosity expansion coefficients"""
            coefs = [ row[Database._CSV_HEADERS[viscosity_key]][0] for viscosity_key in viscosity_keys ]
            ranges = [ row[Database._CSV_HEADERS[viscosity_key]][1] for viscosity_key in viscosity_keys ]
            all_coefs_are_None = [ all([val == None for val in coef_array]) for coef_array in coefs ]
            range_is_None = [ range == None for range in ranges ]

            temperature_range_is_in_wrong_expansion = (all_coefs_are_None[0] and range_is_None[1]) \
                                                    or (all_coefs_are_None[1] and range_is_None[0])

            if all(all_coefs_are_None) and all(range_is_None):
                # No viscosity data is given for this salt, so just remove the second viscosity key and set the other to None
                row.pop(Database._CSV_HEADERS[viscosity_keys[1]])
                row[Database._CSV_HEADERS[viscosity_keys[0]]] = None
            else:
                if temperature_range_is_in_wrong_expansion:
                    # Then switch the temperature range to the expansion with a nonempty coefficient array
                    ranges[0], ranges[1] = ranges[1], ranges[0]
                # Now pop off the empty viscosity expansion
                for viscosity_key, all_coefs_None, coef_array, range in zip(viscosity_keys, all_coefs_are_None, coefs, ranges):
                    if all_coefs_None:
                        row.pop(Database._CSV_HEADERS[viscosity_key])
                    else:
                        row[Database._CSV_HEADERS[viscosity_key]] = [ [0 if coef is None else coef for coef in coef_array], \
                                                                         range]

        def _parse_temp_range(self, col):
            """A function for parsing a temperature range into a tuple or a None"""
            if col == '----' or col == '':
                append_val = None
            else:
                if "-" in col and col.split('-')[0] != '': # The second condition checks for a minus sign
                    if col.split('-')[1] == '': # Upper value not given
                        lower_bound = float(col.split('-')[0])
                        append_val = (lower_bound, lower_bound + 500) # Assume a 500 K temperature range
                    elif col.split('-')[0][-1] == 'E': # A number in scientific notation
                        append_val = float(col)
                    else:
                        append_val = tuple( float(subval) for subval in col.split('-') )
                else:
                    try:
                        append_val = float(col) # An uncertainty, in %
                    except:
                        append_val = None
            return append_val

    def __init__(self, mstdb_tp_path: Path, mstdb_tp_rk_path: Path) -> None:
        """Initializer for the Database class which stores a given version of the MSTDB TP and parses results from it"""
        
        self._parser = self._Parsers(self)
        self.data = self._parser._parse_mstdb_tp(mstdb_tp_path)
        self.rk = self._parser._parse_mstdb_tp_rk(mstdb_tp_rk_path)

    def get_tp(self, thermophysical_property: str, composition_dict: dict) -> Callable[[float], float]:
        """A function for evaluating thermophysical properties for a given salt composition
        using the database
        
        Parameters:
        -----------
            thermophysical_property: Can take the values: {UNIQUE_TP_NAMES}
            composition_dict: A composition dict corresponding to an arbitrary molten salt (with an arbitrary composition)
                              whose endmembers are in the database
        
        Returns:
        --------
            A function that takes the temperature in K and returns the thermophysical property of interest in SI units. For
            clarity, the units are given in the .units attribute of the function
        """

        # ----------------------------------------------------------
        # First perform checks on the input to make sure it's valid
        # ----------------------------------------------------------

        # Verify that the requested thermophysical property is valid
        assert thermophysical_property in Database._UNIQUE_TP_NAMES, (f"The thermophysical property {thermophysical_property} "
                                                                      f"is not one of the valid options: {Database._UNIQUE_TP_NAMES}")

        # Check mole fractions of endmembers add to 1
        assert sum(composition_dict.values()) ==1, ("The mole fractions of the endmembers do not sum to 1, this is an "
                                                    "invalid salt mixture")
        # Make sure endmembers are in database
        database_endmembers = set(key for frozen_dict_key in self.data.keys() for key in frozen_dict_key.keys())
        endmembers_data_dict = {endmember: False for endmember in composition_dict.keys()}
        for endmember in composition_dict.keys():
            assert endmember in database_endmembers, (f"The endmember {endmember} is not contiained in the given database "
                                                       "either update the database path or provide a different salt")
            
            if self.data[frozendict({endmember: 1.0})][thermophysical_property] is not None:
                endmembers_data_dict[endmember] = True
            else:
                # Make sure the user is aware that certain endmembers lack data if they do
                warnings.warn((f"The endmember: {endmember} does not have data (as a pure component) for the "
                               f"requested thermophysical property: {thermophysical_property}, this component "
                               "will be excluded from the property estimation."), UserWarning)

        # Verify that, for at least one of the endmembers, the requested thermophysical property is actually evaluated in the database
        # as a pure compound (not including relevant binary subsystems for now)
        number_of_endmembers_with_data = sum(endmembers_data_dict.values())
        assert number_of_endmembers_with_data != 0, (f"None of the endmembers have any data for the requested thermophysical"
                                                     f"property {thermophysical_property}.")  

        # Check if the database has a binary/ternary composition close to the one given
        input_compositions = np.array([ composition for composition in composition_dict.values() if composition != 0])
        input_composition_keys = [ key for key in composition_dict.keys() if composition_dict[key] != 0]
        keys_with_correct_endmembers = [key for key in self.data.keys() if list(key.keys()) == input_composition_keys ]
        for key_dict in keys_with_correct_endmembers:
            key_compositions = np.array([key_dict[key] for key in input_composition_keys])
            thermophysical_property_exists_in_db = self.data[key_dict][thermophysical_property] is not None
            if np.all(np.isclose(key_compositions, input_compositions, atol=Database.COMPOSITION_TOLERANCE)) and thermophysical_property_exists_in_db:
                return self.data[key_dict][thermophysical_property]

        # ---------------------------------------------------------------------------
        # Now perform ideal property estimation (this will be refined with RK later)
        # ---------------------------------------------------------------------------

        # First, get the thermophysical property functions for the endmembers that have data
        endmembers_with_data = [ endmember for endmember, has_data in endmembers_data_dict.items() if has_data ]
        tp_of_endmembers = []
        mw_endmembers = []
        for endmember in endmembers_with_data:
            database_key = frozendict({endmember: 1.0})
            tp_func = self.data[database_key][thermophysical_property]
            tp_of_endmembers.append(tp_func)

            # Now get molecular weights of endmembers
            mw_endmembers.append(self.data[database_key]['molecular_weight'])

        # Now, rescale compositions of endmembers with data so that they sum to one (essentially assuming the endmembers without
        # data don't exist)
        compositions_with_data = [ composition for endmember, composition in composition_dict.items() \
                                         if endmember in endmembers_with_data ]
        composition_endmembers = [ composition / sum(compositions_with_data) for composition in  compositions_with_data]


        ideal_property = Database.IDEAL_FUNCTIONS[thermophysical_property](tp_of_endmembers, composition_endmembers) \
                          if thermophysical_property != 'density' else \
                          Database.IDEAL_FUNCTIONS[thermophysical_property](tp_of_endmembers, composition_endmembers, mw_endmembers)

        if thermophysical_property == 'density':
            # --------------------------------------------------
            # Perform a redlich kister expansion on the density
            # --------------------------------------------------

            # Get all relevant binary subsystems for which there are rk expansion coefficients
            all_binary_sybsystems = list(combinations(composition_dict, 2)) # Gives a list of pairs of keys, want a list of composition dicts
            
            # Now check which binary subsystems have rk data
            def rk_data_is_available(binary_subsystem):
                database_key = frozenset({*binary_subsystem})
                return database_key in self.rk.keys()
            
            all_binary_sybsystems = [binary_subsystem for binary_subsystem in all_binary_sybsystems \
                                     if rk_data_is_available(binary_subsystem)]

            # Now convert to a dictionary with compositions
            all_binary_sybsystems = [ {endmember: composition_dict[endmember] for endmember in combination } \
                                     for combination in all_binary_sybsystems ]
            
            for binary_subsystem in all_binary_sybsystems:
                # Calculate the contribution to the excess density from each binary subsystem
                excess_density_contribution = self._excess_density_contribution(binary_subsystem)
                ideal_property += excess_density_contribution

        return ideal_property
        

    def _excess_density_contribution(self, binary_subsystem: dict):
        """This function returns a function expressing the contribution of a given binary subsystem to the
        excess density term as a function of temperature"""
        x_1, x_2 = binary_subsystem.values()
        component_1, component_2 = binary_subsystem.keys()
        database_key = frozenset({component_1, component_2})
        coef_array, applicable_range = self.rk[database_key]
        tmin, tmax = applicable_range
        def excess_density_contribution(T: float):
            running_sum = 0
            for index, coefs in enumerate(coef_array):
                A, B = coefs
                L = A + B*T
                running_sum += L*(x_1-x_2)**index
            return x_1*x_2*running_sum*ureg('g/cm^3').to('kg/m^3').magnitude
        # Now set temperature range
        thermo_func = ArbitraryThermoFunction(excess_density_contribution, tmin, tmax, units='kg/m^3')

        return thermo_func
