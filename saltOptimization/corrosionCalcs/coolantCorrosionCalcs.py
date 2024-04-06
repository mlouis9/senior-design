import thermoTools
import thermoToolsAdditions as tta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from thermo.mixture import Mixture
from paths import THERMOCHIMICA, THERMOCHIMICA_FLUORIDE_DATA

# Used for storing outputs in the correct place regardless of where this script is called
script_dir = Path(__file__).parent

#=====================================================================================================
# Path specifications

dataFile = THERMOCHIMICA_FLUORIDE_DATA
scriptName = script_dir / "runThermochimica.ti"

thermochimicaPath = THERMOCHIMICA
outputPath = script_dir / "outputs"

#=====================================================================================================
# Composition parameters

surrogates = {'Am': 'U', 'Cm': 'Nd', 'C': 'Fe', 'Si': 'Fe', 'Mn': 'Fe', 'Np': 'U', 'Mo': 'Fe'}
elements = ['U', 'Pu', 'K', 'Ni', 'Ni', 'Cr', 'Fe', 'Nd', 'Zr', 'Na']

# Density of fuel and clad (kg/m3)
rho_clad = 7870
def rho_coolant(T):
    return 1000*(4.1690 - 9.014E-04*T)

ros = [2E-06, 5E-06] 

fuel_composition = {
    'Pu': 0.289756607,
     'U': 0.2,
    'Am': 0.027444106,
    'Np': 0.023517078,
    'Cm': 0.001890731,
     'K': 0.44999955,
    'Cl': 2.1
}

# Get average molecular weight of the fuel
M_fuel = 0
for element, mole_frac in fuel_composition.items():
    M_fuel += Mixture([element]).MW*mole_frac

stainless_steel_316 = {
    'Cr': 0.181,
    'Ni': 0.113,
    'Mo': 0.014,
    'Mn': 0.02,
    'Si': 0.015,
     'C': 0.0037,
    'Fe': 0.652
}