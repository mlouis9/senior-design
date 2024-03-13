from pathlib import Path
from thermophysicalProperties import Database
import matplotlib.pyplot as plt
import numpy as np
from frozendict import frozendict
import pint
ureg = pint.UnitRegistry(auto_reduce_dimensions=True)
from uncertainties import ufloat

script_dir = Path(__file__).resolve().parent

mstdb_tp_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
mstdb_tp_rk_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties_2.1.0_RK.csv')

db = Database(mstdb_tp_path, mstdb_tp_rk_path)

# ----------------------
# Chloride Test Cases
# ----------------------

# Test case 1: KCl-NaCl

mole_frac_NaCl = np.linspace(0,1,10)
temps = [943, 1033, 1123, 1213, 1303]
for temp in temps:
    test_salts = [ {'KCl': 1-x, 'NaCl': x} for x in mole_frac_NaCl]
    densities = [ db.get_tp('density', test_salt)(temp).nominal_value for test_salt in test_salts ]
    plt.plot(mole_frac_NaCl, densities, label=f"{temp} K")
plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Mole Fraction NaCl')
plt.ylabel('Density kg/m$^3$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'KClNaClCalculated.png' ), dpi=500)
plt.clf()

# Test case 2: KCl-NaCl-UCl3

mole_frac_KCl = np.linspace(0,0.9,10)
temps = [900,985, 1070, 1155, 1240]
for temp in temps:
    test_salts = [ {'KCl': x, 'NaCl': (1-x)*0.25, 'UCl3': (1-x)*0.75} for x in mole_frac_KCl]
    densities = [ db.get_tp('density', test_salt)(temp).nominal_value for test_salt in test_salts ]
    plt.plot(mole_frac_KCl, densities, label=f"{temp} K")
plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Mole Fraction KCl')
plt.ylabel('Density kg/m$^3$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'KClNaClUCl3Calculated.png' ), dpi=500)
plt.clf()

# ----------------------
# Fluoride Test Cases
# ----------------------

# Test case 1: NaF-LiF-ZrF4

mole_frac_ZrF4 = np.linspace(0,1.0,10)
temps = [900,985, 1070, 1155, 1240]
for temp in temps:
    test_salts = [ {'ZrF4': x, 'NaF': (1-x)*0.4, 'LiF': (1-x)*0.6} for x in mole_frac_ZrF4]
    densities = [ db.get_tp('density', test_salt)(temp).nominal_value for test_salt in test_salts ]
    plt.plot(mole_frac_ZrF4, densities, label=f"{temp} K")
plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Mole Fraction ZrF$_4$')
plt.ylabel('Density kg/m$^3$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'NaFLiFZrF4Calculated.png' ), dpi=500)
plt.clf()

# -------------------------
# Heat Capacity Test Cases
# -------------------------

mole_frac_LiF = np.linspace(0,1.0,10)
temp = 1200 # Assume this temp, not clearly stated in the ref
other_endmembers = ['KF', 'NaF']
experiemental_datapoints = {
    'KF': [(0.75, ufloat(69.1, 4.6)), (0.5, ufloat(77.1, 5.8)), (0.25, ufloat(75, 5.5))], 
    'NaF': [(0.8, ufloat(65, 5)), (0.6, ufloat(67, 5)), (0.4, ufloat(70, 6)), (0.2, ufloat(69, 6.5))]
}

for other_endmember in other_endmembers:
    test_salts = [ {'LiF': x, other_endmember: (1-x)} for x in mole_frac_LiF]

    # Need MW of each endmember to convert to heat capacity in a per mole basis
    mw_LiF = db.data[frozendict({'LiF': 1.0})]['molecular_weight']
    mw_other = db.data[frozendict({other_endmember: 1.0})]['molecular_weight']
    mw_test_salts = [ (x*mw_LiF + (1-x)*mw_other)*ureg('g/mol').to('kg/mol').magnitude for x in mole_frac_LiF]

    # Now get the heat capacities and uncertainties
    heat_capacities = [ db.get_tp('liquid_heat_capacity', test_salt)(temp).nominal_value*mw_test_salt \
                       for test_salt, mw_test_salt in zip(test_salts, mw_test_salts) ]
    heat_capacities_unc = [ db.get_tp('liquid_heat_capacity', test_salt)(temp).std_dev*mw_test_salt \
                            for test_salt, mw_test_salt in zip(test_salts, mw_test_salts) ]
    plt.errorbar(mole_frac_LiF, heat_capacities, yerr=heat_capacities_unc, capsize=2.5, elinewidth=1.5, \
                 label=f"LiF-{other_endmember}", markersize=6)
    experimental_unc = [datapoint[1].std_dev for datapoint in experiemental_datapoints[other_endmember]]
    experimental_value = [datapoint[1].nominal_value for datapoint in experiemental_datapoints[other_endmember]]
    xvals = [datapoint[0] for datapoint in experiemental_datapoints[other_endmember]]
    plt.errorbar(xvals, experimental_value, yerr=experimental_unc, fmt='.', capsize=2.5, elinewidth=1.5,\
                 label=f"LiF-{other_endmember} Experimental", markersize=6)

plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Mole Fraction LiF')
plt.ylabel('Heat Capacity $\\frac{\\text{J}}{\\text{K}\\cdot \\text{mol}}$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'LiFCsFetcHeatCapacityCalculated.png' ), dpi=500)
plt.clf()

# --------------------------------
# Thermal Conductivity Test Cases
# --------------------------------

# KCl-MgCl2 (0.68-0.32)
test_salt = {'KCl': 0.68, 'MgCl2': 0.32}
temps = np.linspace(450 + 273.15,850 + 273.15, 10)
thermal_conductivities = [db.get_tp('thermal_conductivity', test_salt)(temp).nominal_value for temp in temps]
thermal_conductivities_unc = [db.get_tp('thermal_conductivity', test_salt)(temp).std_dev for temp in temps]
plt.errorbar(temps - 273.15, thermal_conductivities, yerr=thermal_conductivities_unc, capsize=2.5)
# plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Temperature ($^\circ \\text{C}$)')
plt.ylabel('Thermal Conductivity $\\frac{\\text{W}}{\\text{m}\\cdot \\text{K}}$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'MgCl2KClThermalConductivity.png' ), dpi=500)
plt.clf()

# ----------------------
# Viscosity Test Cases
# ----------------------

test_salts = [ 
    {'LiF': 1/3,'NaF': 1/3, 'KF': 1/3} ,
    {'LiF': 0.66, 'BeF2': 0.34},
    {'LiF': 0.73, 'BeF2': 0.27}
]
labels = ['FLiNaK', 'LiF-BeF2 (66-34 mol%)', 'LiF-BeF2 (73-27 mol%)']
temps = np.linspace(600 + 273.15, 700 + 273.14, 10)
for index, test_salt in enumerate(test_salts):
    viscosity = [db.get_tp('viscosity', test_salt)(temp).nominal_value*1E+03 \
                 for temp in temps]
    viscosity_unc = [db.get_tp('viscosity', test_salt)(temp).std_dev*1E+03 \
                     for temp in temps]
    plt.errorbar(temps - 273.15, viscosity, yerr=viscosity_unc, capsize=2.5, elinewidth=1.5, label=labels[index])
plt.legend(loc='upper left', bbox_to_anchor=(0.9,1))
plt.xlabel('Temperature ($^\circ \\text{C}$)')
plt.ylabel('Dynamic Viscosity $\\text{mPa}\cdot \\text{s}$')
plt.grid()
plt.tight_layout()
plt.savefig(str( script_dir / 'plots' / 'FLiNaKViscosity.png' ), dpi=500)
plt.clf()