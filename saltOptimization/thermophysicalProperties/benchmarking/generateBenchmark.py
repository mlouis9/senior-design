from pathlib import Path
from thermophysicalProperties import Database
import matplotlib.pyplot as plt
import numpy as np

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