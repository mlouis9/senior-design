

# Notes
This is true in general with Jupyter notebooks, but it's a good idea to try restarting the jupyter kernel (which deletes any saved global variables) if you're met with any errors when trying to run the notebook for example
```
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Cell In[7], line 26
     23 mstdb_tp_path = Path('../../mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
     24 mstdb_tp_rk_path = Path('../../mstdb-tp/Molten_Salt_Thermophysical_Properties_RK.csv')
---> 26 db = Database(mstdb_tp_path, mstdb_tp_rk_path)

File ~/PythonProjects/seniorDesign/modules/thermophysicalProperties.py:601, in Database.__init__(self, mstdb_tp_path, mstdb_tp_rk_path)
    598 """Initializer for the Database class which stores a given version of the MSTDB TP and parses results from it"""
    600 self._parser = self._Parsers(self)
--> 601 self.data = self._parser._parse_mstdb_tp(mstdb_tp_path)
    602 self.rk = self._parser._parse_mstdb_tp_rk(mstdb_tp_rk_path)

File ~/PythonProjects/seniorDesign/modules/thermophysicalProperties.py:379, in Database._Parsers._parse_mstdb_tp(self, mstdb_tp_path)
    377 for tp_key in [ key for key in Database._TP_FUNCTIONS.keys() if key in parsed_row.keys() ]:
    378     if parsed_row[tp_key] is not None: # Can be none if tp_key is the wrong viscosity functionalization
--> 379         parsed_row[tp_key] = ThermoFunction(parsed_row, tp_key)
    380     if 'viscosity' in tp_key:
    381         parsed_row['viscosity'] = parsed_row.pop(tp_key)

File ~/PythonProjects/seniorDesign/modules/thermophysicalProperties.py:153, in ThermoFunction.__init__(self, salt, key)
    150 else:
    151     key_of_thermo_func = self.key
--> 153 if isinstance(salt[key_of_thermo_func], ArbitraryThermoFunction):
    154     # If salt[key_of_thermo_func] is already a ThermoFunction, no need to parse
    155     # they min/max temp, etc. just take it from the ThermoFunction attributes
    156     func = salt[key_of_thermo_func].func
    157     min_temp = salt[key_of_thermo_func].min_temp

KeyError: 'viscosity'
```
which causes trouble because the database instance is already partially parsed (as it was parsed previously) and so it doesn't have the expected structure and the parsing script fails.