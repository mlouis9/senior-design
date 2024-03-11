import csv
from pathlib import Path
from typing import Dict, Any
from io import StringIO
from frozendict import frozendict

""" This is a module for reading and calculating thermophysical properties from the MSTDB using ideal estimations (e.g. 
additivity of molar volumes) and the RK expansion for estimating the effects nonideal mixing.
"""

class Database:
    def __init__(self, mstdb_tp_path: Path, CSV_HEADERS: dict=None, CSV_VALUES: dict=None) -> None:
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
        self.data = self._parse_mstdb_tp(mstdb_tp_path)

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
    
    def _preprocess_subheader(self, reader, csvfile):
        """Function for correctly processing the subheaders in MSTDB-TP corresponding to the coefficients of the temperature
        dependence expansion of TP's"""
        # First remove the first row of reader, it corresponds to the subheader
        reader.pop(0)
        # reader = [{key: value for key, value in row.items() if key != ''} for row in reader]

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
                    print(reader[row_index-3][key_of_last_nonempty_col][0])
                    reader[row_index-3][key_of_last_nonempty_col][0].append(append_val)
                else:
                    if index_of_last_empty_col == col_index - 1:
                        # Add temp range (a string) to the end of the list
                        if col == '----' or col == '':
                            append_val = None
                        else:
                            append_val = col
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
        if all( [ row[self._CSV_HEADERS[viscosity_key]] == None for viscosity_key in viscosity_keys ] ):
            # Neither viscosity key is given, just for completeness take the first key (even though it's None)
            row.pop(self._CSV_HEADERS[viscosity_keys[1]])
        else:
            # Git rid of the other key
            for viscosity_key in viscosity_keys:
                if row[self._CSV_HEADERS[viscosity_key]] == None:
                    row.pop(self._CSV_HEADERS[viscosity_key])
        # Now convert all keys to the convenient names
        def get_first_key_with_value(value_to_match, my_dict):
            return next( ( key for key, value in my_dict.items() if value ==  value_to_match), None )
        # parsed_row = {get_first_key_with_value(key, self._CSV_HEADERS): value for key, value in parsed_row.items() }
        print(parsed_row)
        return parsed_row

# Example usage:
mstdb_tp_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
db = Database(mstdb_tp_path)
# print(db.data)  # This will print the parsed data as a dictionary with frozendict keys