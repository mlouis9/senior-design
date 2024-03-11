import csv
from pathlib import Path
from typing import Dict, Any
from io import StringIO
from frozendict import frozendict

""" This is a module for reading and calculating thermophysical properties from the MSTDB using ideal estimations (e.g. 
additivity of molar volumes) and the RK expansion for estimating the effects nonideal mixing.
"""

class Database:
    def __init__(self, mstdb_tp_path: Path, mstdb_tp_rk_path: Path, CSV_HEADERS: dict=None, CSV_VALUES: dict=None, \
                 CSV_HEADERS_RK: dict=None) -> None:
        """Initializer for the Database class which stores a given version of the MSTDB TP and parses results from it"""
        if CSV_HEADERS == None:
            # Configuration for CSV headers
            self._CSV_HEADERS = {
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
        if CSV_VALUES == None:
            self._CSV_VALUES = {
                'pure_salt': 'Pure Salt'
            }
        if CSV_HEADERS_RK == None:
            # Headers for the RK expansion coefficient file
            self._CSV_HEADERS_RK = {
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
        self.data = self._parse_mstdb_tp(mstdb_tp_path)
        self.rk = self._parse_mstdb_tp_rk(mstdb_tp_rk_path)

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
                data.update({composition_dict: parsed_row})
        return data
    
    def _parse_mstdb_tp_rk(self, mstdb_tp_rk_path):
        """Function for parsing the .csv file containing the RK coefficients for the density expansion. Since this
        csv does not have a subheader, its processing should be much easier"""
        rk = dict()
        print(mstdb_tp_rk_path)
        with mstdb_tp_rk_path.open('r', errors='replace') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')
            for row in reader:
                salt_pair_key = frozenset({row[self._CSV_HEADERS_RK['component1']],\
                                           row[self._CSV_HEADERS_RK['component2']]})
                parsed_row = self._parse_row_rk(row)
                rk.update( { salt_pair_key: parsed_row } )
                
        return rk
    
    def _parse_row_rk(self, row):
        """Function for parsing rows of the MSTDB-TP RK csv"""
        parsed_row = [[],]
        keys_by_coefficient_order = [ ( self._CSV_HEADERS_RK[f'a{i}'], self._CSV_HEADERS_RK[f'b{i}'] ) for i in range(1,4) ]

        # First add coefficients for each order
        for keys_at_order in keys_by_coefficient_order:
            parsed_row[0].append([float(row[keys_at_order[0]]), float(row[keys_at_order[1]])])

        # Then add temperature range
        parsed_row.append((float(row[self._CSV_HEADERS_RK['tmin']]), float(row[self._CSV_HEADERS_RK['tmax']])))

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
        for row_index, row in enumerate(csvfile):
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
                        # Add temp range to the end of the list
                        append_val = self._parse_temp_range(col)
                            
                        reader[row_index -3][key_of_last_nonempty_col][1] = append_val

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
        if row[self._CSV_HEADERS['composition']] == self._CSV_VALUES['pure_salt']:
            # The composition of the salt is just 100% of the listed formula for pure salts
            return frozendict({row[self._CSV_HEADERS['formula']]: 1.0})
        else:
            composition = row[self._CSV_HEADERS['composition']].split('-')
            formula = row[self._CSV_HEADERS['formula']]
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
            if key in [self._CSV_HEADERS[col_header] for col_header in cols_to_skip  ]:  # Skip the composition column as it is handled separately
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
        parsed_row = {get_first_key_with_value(key, self._CSV_HEADERS): value for key, value in parsed_row.items() }
        return parsed_row
    
    def _parse_viscosity_coefs(self, viscosity_keys, row):
        """This function is meant for processing viscosity expansion coefficients"""
        coefs = [ row[self._CSV_HEADERS[viscosity_key]][0] for viscosity_key in viscosity_keys ]
        ranges = [ row[self._CSV_HEADERS[viscosity_key]][1] for viscosity_key in viscosity_keys ]
        all_coefs_are_None = [ all([val == None for val in coef_array]) for coef_array in coefs ]
        range_is_None = [ range == None for range in ranges ]

        temperature_range_is_in_wrong_expansion = (all_coefs_are_None[0] and range_is_None[1]) \
                                                or (all_coefs_are_None[1] and range_is_None[0])

        if all(all_coefs_are_None) and all(range_is_None):
            # No viscosity data is given for this salt, so just remove the viscosity keys
            for viscosity_key in viscosity_keys:
                row.pop(self._CSV_HEADERS[viscosity_key])
        else:
            if temperature_range_is_in_wrong_expansion:
                # Then switch the temperature range to the expansion with a nonempty coefficient array
                ranges[0], ranges[1] = ranges[1], ranges[0]
            # Now pop off the empty viscosity expansion
            for viscosity_key, all_coefs_None, coef_array, range in zip(viscosity_keys, all_coefs_are_None, coefs, ranges):
                if all_coefs_None:
                    row.pop(self._CSV_HEADERS[viscosity_key])
                else:
                    row[self._CSV_HEADERS[viscosity_key]] = [coef_array, range]

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
                append_val = col
        return append_val

# Example usage:
mstdb_tp_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
mstdb_tp_rk_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties_2.1.0_RK.csv')

db = Database(mstdb_tp_path, mstdb_tp_rk_path)
example_salt = frozendict({'AlCl3': 1.0})
# print(db.data[example_salt])  # This will print the parsed data as a dictionary with frozendict keys
print(db.rk)