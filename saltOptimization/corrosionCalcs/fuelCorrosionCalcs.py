import thermoTools
import thermoToolsAdditions as tta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from thermo.mixture import Mixture
import os, sys

# Get the script name from sys.argv
script_name = sys.argv[0]

# Get the absolute path of the script
script_dir = Path(os.path.abspath(script_name)).resolve().parent

#=====================================================================================================
# Path specifications

dataFile = script_dir / "../../thermochimica/data/MSTDB-TC_V3.0_Chlorides_No_Functions_8-2.dat"
scriptName = script_dir / "runThermochimica.ti"

thermochimicaPath = script_dir / "../../thermochimica"
outputPath = script_dir / "outputs"

#=====================================================================================================
# Composition parameters

surrogates = {'Am': 'U', 'Cm': 'Nd', 'C': 'Fe', 'Si': 'Fe', 'Mn': 'Fe', 'Np': 'U', 'Mo': 'Fe'}
elements = ['U', 'Pu', 'K', 'Cl', 'Ni', 'Cr', 'Fe', 'Nd', 'Zr']

# Density of fuel and clad (kg/m3)
rho_clad = 7870
def rho_fuel(T):
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

# Replace elements in fuel and structural material with surrogates
for element, surrogate in surrogates.items():
    # Replace element by its surrogate
    for composition_dict in [fuel_composition, stainless_steel_316]:
        if element in composition_dict:
            # Check if surrogate key already exists, if so, update it
            if surrogate in composition_dict.keys():
                composition_dict[surrogate] += composition_dict[element]
            else:
                composition_dict.update({surrogate: composition_dict[element]})
            composition_dict.pop(element) # Remove old element
#=====================================================================================================

tstart = 900
tend = 1500
ntstep = 10
pstart = 1
pend = 1
npstep = 1

#=====================================================================================================

zr_fractions = [0.0, 1.0]

# Initialize dictionaries that will store the results of the calculations
corroded_element_moles = dict()
n_corr = dict()
M_corr = dict()
colors = {'No Zr': 'tab:orange', 'Zr': 'tab:blue'}
for zr_fraction in zr_fractions:
    zr_composition = {'Zr': zr_fraction}

    # Now add the dictionaries together
    total_composition = Counter(fuel_composition) + Counter(stainless_steel_316)

    # Now convert back to dictionary
    total_composition = dict(total_composition)
    total_composition.update({'Zr': zr_fraction})

    masses = [ total_composition[element] for element in elements ]

    # Write an input script
    thermoTools.WriteInputScript(scriptName, str(dataFile), elements, tstart, tend, ntstep, pstart, pend, npstep, masses, \
                                    fuzzyStoichiometry=True, gibbsMinCheck=True)
    # Run an input script
    outputName = f"corrZr{zr_fraction:1.0f}.json"

    thermoTools.RunInputScript(scriptName, jsonName=str(outputPath / outputName), thermochimica_path=str(thermochimicaPath),\
                                noOutput=True)

    # Now read the output using tta and plot
    calc = tta.thermoOut(str(script_dir / f'outputs/{outputName}'))

    # Make sure that no surrogate elements are also present in the fuel (otherwise we have to add some new capabilities)
    for element in stainless_steel_316.keys():
        assert element not in fuel_composition.keys()

    if zr_fraction == 0:
        label = "No Zr" 
    else:
        label = f"Zr"

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
    for index, ro in enumerate(ros):
        corrosion_depth =1/2*M_corr[label]/M_fuel*rho_fuel(temperatures)/rho_clad*ro*n_corr[label]
        if index == 1:
            linestyle = 'dashed'
        else:
            linestyle = 'solid'

        plt.plot(temperatures, corrosion_depth/1E-09, label=f"{label} ($r_o$={ro*1E+06}mm)", marker='.', linestyle=linestyle, \
                 color = colors[label])

# Now finish off plot and save
plt.yscale('log')
plt.ylabel('Total corrosion depth ($\mu$m)')
plt.xlabel('Temperature (K)')
plt.grid()
plt.legend()
plt.savefig(str(script_dir / "plots/effectOfZr.png"), dpi=500)
plt.clf()
            
# Plot corrosion composition
for label, corroded_elements in corroded_element_moles.items():
    for element, mole_corroded in corroded_elements.items():
        plt.plot(temperatures, mole_corroded/n_corr[label]*100, label=f"{element} ({label})", marker='.')
    plt.yscale('log')
    plt.ylabel('Fraction of element dissolved in fuel salt (%)')
    plt.xlabel('Temperature (K)')
    plt.grid()
    plt.legend()
    plt.savefig(str(script_dir / f"plots/corrodedElements{label}.png"), dpi=500)
    plt.clf()

# Now plot mol CrCl2/mol fuel
for label, corroded_elements in corroded_element_moles.items():
    n_cr = corroded_elements['Cr']
    plt.plot(temperatures, n_cr, label=label, marker='.', color=colors[label])
plt.yscale('log')
plt.ylabel('moles of Cr corroded / mole of fuel salt')
plt.xlabel('Temperature (K)')
plt.grid()
plt.legend()
plt.savefig(str(script_dir / 'plots/chromiumVsZr.png'), dpi=500)
plt.clf()