import thermoTools
import modules.thermoToolsAdditions as tta
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# File IO input parameters
thermochimica_path = Path("/home/mlouis9/thermochimica")
output_path = Path("../../PythonProjects/fuelOptimization/outputs")
output_name = 'output.json'
data_file = Path("../data/MSTDB-TC_V3.0_Chlorides_No_Functions_8-2.dat")
input_file_name = "runThermochimica.ti"


# Physical input parameters (same for all pseudo binary calcaultions)
press = 100
tunit = 'K'
punit = 'atm'
munit = 'moles'

xlo = 0.0
xhi = 1.0
nxstep = 120
tlo = 500
thi = 1200
ntstep = 120

# Define cases, i.e. all of the pseudo binary systems shown in Zhang
cases = [(['K', 'Cl', 'Pu'], {'K Cl': 1.0}, {'Pu Cl_3': 1.0}), \
         (['Li', 'Cl', 'Pu'], {'Li Cl': 1.0}, {'Pu Cl_3': 1.0}), \
         (['K', 'Cl', 'Li'], {'Li Cl': 1.0}, {'K Cl': 1.0}), \
         (['Pu', 'Cl', 'K', 'Li'], {'K Cl': 0.414, 'Li Cl': 0.586}, {'Pu Cl_3': 1.0})] # Solubility calculation

case_names = ["KCl-PuCl3", "LiCl-PuCl3", "LiCl-KCl", "KClLiCl-PuCl3"]

for case_index, case in enumerate(cases):
    elements_used = case[0]
    left_endmember_composition = case[1]
    right_endmember_composition = case[2]

    # Run thermochimica calculation
    calc = tta.pseudo_binary_calculation(thermochimica_path, output_path, output_name, data_file, xlo, xhi, nxstep, tlo, thi, ntstep, elements_used,\
                                left_endmember_composition, right_endmember_composition)

    # Now load output.json and plot

    # Plot region diagram
    loaded_output = tta.thermoOut(Path('../outputs/output.json'))
    diagram = tta.pseudoBinaryDiagram(left_endmember_composition, right_endmember_composition, Path('../outputs/output.json'), \
                                    plot_everything=case_index == 3, ntstep=ntstep, nxstep=nxstep)
    diagram.plot_phase_regions(plot_marker='.', plot_mode='region')
    diagram.plot.fig.savefig(f"plots/{case_names[case_index]}Region", bbox_inches='tight')

    # Plot boundary diagram
    diagram.plot_phase_regions(plot_marker='.-', plot_mode='boundary')
    diagram.plot.fig.savefig(f"plots/{case_names[case_index]}Boundary", bbox_inches='tight')

    # If solubility calc, calculate solubility as a function of temperature and plot
    if case_index == 3:
        plt.clf()
        temps = np.linspace(tlo, thi, ntstep)
        solubility = diagram.calculate_solubility(frozenset({'MSCL', 'PuCl3_P63/M_No.176(s)'}), temps)
        nonzero_indices = np.nonzero(solubility)
        plt.plot(temps[nonzero_indices], solubility[nonzero_indices], label='Thermochimica')

        ref_solubility = np.array([0.382, 0.415, 0.460])
        ref_temps = np.array([723, 773, 823])
        plt.scatter(ref_temps, ref_solubility, label='Reference')
        plt.grid()
        plt.ylim(0.3, 0.6)
        plt.xlim(700, 850)
        plt.xlabel('T (K)')
        plt.ylabel('$S_{\\text{PuCl}_3}\\text{ (mole fraction PuCl}_3\\text{)}$')
        plt.legend()
        plt.savefig('plots/solubilityComparison')
