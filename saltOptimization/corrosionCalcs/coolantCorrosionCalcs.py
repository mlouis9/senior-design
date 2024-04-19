import thermoTools
import thermoToolsAdditions as tta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from thermo.mixture import Mixture
from paths import THERMOCHIMICA, THERMOCHIMICA_FLUORIDE_DATA, MSTDB_TP_DATA, MSTDB_TP_RK_DATA
from frozendict import frozendict
import warnings

warnings.filterwarnings('ignore')

# Get database for thermophysical properties
from thermophysicalProperties import Database

db = Database(MSTDB_TP_DATA, MSTDB_TP_RK_DATA)

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

surrogates = {'C': 'Fe', 'Si': 'Fe', 'Mn': 'Fe', 'Mo': 'Fe'}
elements = ['F', 'Na', 'Zr', 'K', 'Li', 'Ni', 'Cr', 'Fe']

# Density of coolant and clad (kg/m3)
rho_clad = 7870

coolant = frozendict({'NaF': 0.132, 'ZrF4': 0.032, 'KF': 0.4109, 'LiF': 0.4251})
rho_coolant = db.get_tp('density', coolant, uncertainty=False)

# Now the composition of the coolant must be defined for thermochimica calculations
coolant_composition = {
    'Na F': 0.132, 
    'Zr F_4': 0.032, 
    'K F': 0.4109, 
    'Li F': 0.4251
}
unique_elements = ['F', 'Na', 'Zr', 'K', 'Li']

# Now parse into mass labels
coolant_composition = np.array(tta.component_fractions_to_element_fractions(coolant_composition, unique_elements))
coolant_composition /= np.sum(coolant_composition)
coolant_composition = {element: fraction for element, fraction in zip(unique_elements, coolant_composition)}

# Get average molecular weight of the coolant
M_coolant = 0
for element, mole_frac in coolant_composition.items():
    M_coolant += Mixture([element]).MW*mole_frac

stainless_steel_316 = {
    'Cr': 0.181,
    'Ni': 0.113,
    'Mo': 0.014,
    'Mn': 0.02,
    'Si': 0.015,
     'C': 0.0037,
    'Fe': 0.652
}

# Replace surrogates for Stainless steel (each of the endmembers for the coolant salt are in the database)
for element, surrogate in surrogates.items():
    # Replace element by its surrogate
    if element in stainless_steel_316:
        # Check if surrogate key already exists, if so, update it
        if surrogate in stainless_steel_316.keys():
            stainless_steel_316[surrogate] += stainless_steel_316[element]
        else:
            stainless_steel_316.update({surrogate: stainless_steel_316[element]})
        stainless_steel_316.pop(element) # Remove old element

#=====================================================================================================
# Thermochimica calculation parameters

tstart = 750
tend = 1000
ntstep = 10
pstart = 1
pend = 1
npstep = 1

#=====================================================================================================

ro = 5.5E-03 # m

# Initialize dictionaries that will store the results of the calculations
corroded_element_moles = dict()
n_corr = dict()
M_corr = dict()

ZrF2_fractions = np.linspace(0, 1, 5) # Note, cannoet exceed ZrF4 Composition in coolant salt
colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(ZrF2_fractions))))

for ZrF2_fraction in ZrF2_fractions:

    # Make coolant composition
    coolant_composition = {
        'Na F': 0.132, 
        'Zr F_4': 0.032*(1-ZrF2_fraction),
        'Zr F_2': 0.032*ZrF2_fraction,
        'K F': 0.4109, 
        'Li F': 0.4251
    }

    unique_elements = ['F', 'Na', 'Zr', 'K', 'Li']

    # Now parse into mass labels
    coolant_composition = np.array(tta.component_fractions_to_element_fractions(coolant_composition, unique_elements))
    coolant_composition /= np.sum(coolant_composition)
    coolant_composition = {element: fraction for element, fraction in zip(unique_elements, coolant_composition)}

    # Get average molecular weight of the coolant
    M_coolant = 0
    for element, mole_frac in coolant_composition.items():
        M_coolant += Mixture([element]).MW*mole_frac

    # Now add together the cooalnt and stainless steel compositions to get the total composition for the thermochimica calc
    total_composition = Counter(coolant_composition) + Counter(stainless_steel_316)

    # Now convert back to dictionary
    total_composition = dict(total_composition)

    masses = [ total_composition[element] for element in elements ]

    # --------------------------
    # Thermochimica Calculation
    # --------------------------
    
    # Write an input script
    thermoTools.WriteInputScript(scriptName, str(dataFile), elements, tstart, tend, ntstep, pstart, pend, npstep, masses, \
                                    fuzzyStoichiometry=True, gibbsMinCheck=True)
    # Run an input script
    outputName = f"corrZrF2-{ZrF2_fraction*100:1.0f}.json"

    thermoTools.RunInputScript(scriptName, jsonName=str(outputPath / outputName), thermochimica_path=str(thermochimicaPath), \
                                    noOutput=True)

    # --------------------
    # Postprocess Results
    # --------------------
    
    # Now read the output using tta and plot
    calc = tta.thermoOut(str(script_dir / f'outputs/{outputName}'))

    label = f"ZrF$_2$/ZrF$_4$ {ZrF2_fraction*100:2.1f}%"

    nstates = len(list(calc.output.keys()))
    M_corr.update({label: np.zeros(nstates)})
    n_corr.update({label: np.zeros(nstates)})

    # Corroded elements that we want to keep track of

    corroded_elements = ['Cr', 'Fe', 'Ni']
    corroded_element_moles.update({label: {element: np.zeros(nstates) for element in corroded_elements}})
    temperatures = np.zeros(nstates)

    # Iterate over states in output
    for state_index, state in enumerate(calc.output.keys()):
        # now iterate over elements in stainless_steel_316 material
        for element, mole_frac in stainless_steel_316.items():
            # Fraction of element in solution
            soln_frac_element = calc.solution_fraction_element[element][state]
            n_soln_element = soln_frac_element*mole_frac

            n_corr[label][state_index] += n_soln_element
            M_corr[label][state_index] += np.multiply(n_soln_element, calc.atomic_weights[element])
        
        # Get number of moles of chromium, for comparison
        for element in corroded_elements:
            corroded_element_moles[label][element][state_index] = calc.solution_fraction_element[element][state]\
                                                            *stainless_steel_316[element]

        # Now get temperature at this state
        temperatures[state_index] = calc.output[state]['temperature']

    # Now divide M_corr by n_corr to get average molecular weight
    M_corr[label] = M_corr[label]/n_corr[label]

    # Now calculate the corrosion depth
    corrosion_depth =1/2*M_corr[label]/M_coolant*rho_coolant(temperatures)/rho_clad*ro*n_corr[label]

    plt.plot(temperatures, corrosion_depth/1E-06, label=f"{label} ($r_o$={ro*1E+03}mm)", marker='.', \
                color = next(colors))

# Now finish off plot and save
plt.yscale('log')
plt.ylabel('Total corrosion depth ($\mu$m)')
plt.xlabel('Temperature (K)')
plt.grid()
plt.legend()
plt.savefig(str(script_dir / "plots/effectOfZrF2Coolant.png"), dpi=500)
plt.clf()