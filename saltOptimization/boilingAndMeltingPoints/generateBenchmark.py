import thermoToolsAdditions as tta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from paths import THERMOCHIMICA, THERMOCHIMICA_FLUORIDE_DATA

script_dir = Path(__file__).parent

# File IO input parameters
tta.Config.set_config(
    thermochimica_path = THERMOCHIMICA,
    output_path = script_dir,
    data_file = THERMOCHIMICA_FLUORIDE_DATA,
)

salt_composition = {'Li F':0.465, 'Na F': 0.115, 'K F': 0.42}
elements_used = ['Li', 'F', 'Na', 'K']

# NOTE using default tlo=0, thi=2500, atmospheric pressure, MSCL as the liquid phase, and gas_ideal as the gaseous phase, these should
# be appropriate for all calculations, but are provided as arguments for flexibility.

T_m, T_b = tta.calculate_melting_and_boiling(salt_composition, elements_used, method='phase diagram')
print(f"Phase diagram method {T_m}, {T_b}")

# Try this using the faster calclist method
T_m, T_b = tta.calculate_melting_and_boiling(salt_composition, elements_used, phase_tolerance=0)
print(f"Fast method {T_m}, {T_b}")

# Using a differrent phase tolerance (the default is 0.9)
T_m, T_b = tta.calculate_melting_and_boiling(salt_composition, elements_used, phase_tolerance=0)
print(f"Fast method with a nonzero phase tolerance {T_m}, {T_b}")

output = tta.thermoOut(default=True)
temperatures = output.temperatures.values()
mole_fraction_liquid = output.get_phase_fraction('MSFL')
mole_fraction_liquid_3 = output.get_phase_fraction('MSFL#3')
mole_fraction_gas = output.get_phase_fraction('gas_ideal')

plt.plot(temperatures, np.array(mole_fraction_liquid) + np.array(mole_fraction_liquid_3), label='MSFL')
# plt.plot(temperatures, mole_fraction_liquid_3, label='MSFL #3')
plt.plot(temperatures, mole_fraction_gas, label='Gas')
plt.xlabel('Temperature (K)')
plt.ylabel('Mole Fraction of Phase')
plt.grid()
plt.legend()
plt.savefig(str(script_dir / 'phaseFractions.png'), dpi=500)