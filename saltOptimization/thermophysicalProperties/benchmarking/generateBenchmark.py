from pathlib import Path
from thermophysicalProperties import Database
import matplotlib.pyplot as plt
import numpy as np

script_dir = Path(__file__).resolve().parent

mstdb_tp_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
mstdb_tp_rk_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties_2.1.0_RK.csv')

db = Database(mstdb_tp_path, mstdb_tp_rk_path)

# ----------------------
# Test case 1: KCl-NaCl
# ----------------------

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