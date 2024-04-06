from pathlib import Path
from thermophysicalProperties import Database
import matplotlib.pyplot as plt
import numpy as np
from frozendict import frozendict
import pint
ureg = pint.UnitRegistry(auto_reduce_dimensions=True)
from uncertainties import ufloat
import os, sys
from numdifftools import Derivative
from scipy.optimize import minimize
from paths import ROOT, THERMOCHIMICA_FLUORIDE_DATA

print(THERMOCHIMICA_FLUORIDE_DATA)
# warnings.filterwarnings("ignore")

# # Get the script name from sys.argv
# script_name = sys.argv[0]

# # Get the absolute path of the script
# script_dir = Path(os.path.abspath(script_name)).resolve().parent

# mstdb_tp_path = Path('../../mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
# mstdb_tp_rk_path = Path('../../mstdb-tp/Molten_Salt_Thermophysical_Properties_RK.csv')

# db = Database(mstdb_tp_path, mstdb_tp_rk_path)

# import thermoToolsAdditions as tta
# from thermoToolsAdditions import convert_to_thermochem_name

# # File IO input parameters
# thermochimica_path = Path("../../thermochimica")
# output_path = Path('../outputs')
# output_name = 'output.json'
# data_file = Path("../../thermochimica/data/MSTDB-TC_V3.0_Chlorides_No_Functions_8-2.dat")
# input_file_name = "runThermochimica.ti"

# # First do the melting/boiling point calculation (this performs a series of thermochimica calculations)
# test_composition = {'Na Cl': 0.3, 'K Cl': 0.05, 'Zr Cl_4': 0.1, 'Pu Cl_3': 0.398, 'Nd Cl_3': 0.002, 'U Cl_3': 0.15} #'K Cl': 0.01, 'Al Cl_3': 0.01, 'Zr Cl_4': 0.28
# elements_used = ['Cl', 'Na', 'K', 'Zr', 'Pu', 'Nd', 'U']
# print(tta.calculate_melting_and_boiling(thermochimica_path, output_path, output_name, data_file, \
#                                              test_composition, elements_used, ntstep=150))
# output = tta.thermoOut(output_path / output_name)
# temperatures = output.temperatures.values()
# mole_fraction_liquid = output.get_phase_fraction('MSCL')
# mole_fraction_liquid_3 = output.get_phase_fraction('MSCL#3')
# mole_fraction_gas = output.get_phase_fraction('gas_ideal')

# plt.plot(temperatures, np.array(mole_fraction_liquid) + np.array(mole_fraction_liquid_3), label='MSCL')
# plt.plot(temperatures, mole_fraction_liquid_3, label='MSCL #3')
# plt.plot(temperatures, mole_fraction_gas, label='Gas')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Mole Fraction of Phase')
# plt.grid()
# plt.legend()
# plt.show()
