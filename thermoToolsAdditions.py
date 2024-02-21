import thermoTools
import json
from pathlib import Path
import numpy as np
import os

"""This is a module containing all of the supporting classes and functions for running the calculations needed for my Senior Deisgn
project. This module extends the built-in thermochimica module `thermoTools`, and contains a function for automatically executing a
solubility calculation, and easily reading output.

Author: Matthew Louis
Email:" matthewlouis31@gmail.com
"""

# ---------------------------------
#           Classes
# ---------------------------------

class thermoOut:
    """A class for storing Thermochimica output, with methods for easily plotting
    and postprocessing"""
    def __init__(self, out_file: Path=None):
        """Initialize a thermoOut object
        
        Parameters
        ----------
            outFile: The absolute path (or path relative to the python script utilizing this) to the
                     Thermochimica output .json file
        
        Returns
        -------
            Nothing, class initializer
        """

        # If an outFile is not provided, do nothing
        if out_file != None:
            # First read .json and save as .output attribute
            with open(str(out_file), 'r') as f:
                self.output = json.load(f)

            # Read the temperatures and pressures represented in this calculation
            self._reinitialize()
        else:
            self.output = {}
            

    def _reinitialize(self, null=False):
        """Uses internal class methods which manipulate the output json data to initialize important derived quantities/data
        in object attributes.
        
        Parameters:
        -----------
            null: A flag indicating whether the reinitialization should take place with or without output (this is for supporting
                  the creation of 'empty' objects, i.e. my_object = thermoOut()), if True, initialization without output, and thus
                  all relevant attributes are assigned empty objects
        
        Returns:
        --------
            Nothing
        """
        # Read the temperatures and pressures represented in this calculation
        if not null:
            # self.temperatures = np.array([ state['temperature'] for state in self.output.values()])
            arr = []
            for key, state in self.output.items():
                try:
                    arr.append(state['temperature'])
                except:
                    print(key)
            self.temperatures = np.array(arr)

            self.stable_phases = self._get_stable_phases()
            self.elements = set( element for state in  list(self.output.values()) \
                                 for element in list(state['elements'].keys()) ) # Get unique elements corresponding to each of the
                                                                                 # states, should be the SAME for all states (although
                                                                                 # different mol fractions are possible). NOTE that
                                                                                 # thermochimica does not add elements with zero mole
                                                                                 # fractions to output.json, though they should still
                                                                                 # be included for our purposes
            self.mole_fraction_element_by_phase = self._get_mole_fraction_element_by_phase()

        else:
            self.temperatures = []
            self.stable_phases = []
            self.mole_fraction_element_by_phase = {}
            self.elements = []

    def add_output(self, out_file):
        """NOTE: This function does NOT check to see if states are unique before adding them to the output, it is possible
        to end up with duplicate states if this method is not used judiciously"""

        # Read the new outfile
        with open(out_file, 'r') as f:
            new_output = json.load(f)

        # Now find the last state in the current output
        if self.output != {}:
            last_state = int(list(self.output.keys())[-1]) # Note the state keys are 1-indexed
        else:
            last_state = 0 # If self.output = {}, list(self.output.keys()) has no -1th element, so the above errors

        for index, state in enumerate(new_output.values()):
            self.output.update({str(last_state + index + 1): state})

        self._reinitialize()


    def _get_stable_phases(self, DEBUG=False):
        # Below this tolerance, set phase fraction = 0
        phase_include_tol = 1e-8

        stable_phases = []
        for state in self.output.values():
            solution_phases = [ (name, phase['moles']) for name, phase in  state['solution phases'].items() \
                               if phase['moles'] > phase_include_tol ]
            pure_condensed_phases = [ (name, phase['moles']) for name, phase in  state['pure condensed phases'].items() \
                                     if phase['moles'] > phase_include_tol ]
            if DEBUG and (len(solution_phases) == 0 or 'MSCL' not in [solution_phase[0] for solution_phase in solution_phases]):
                print([solution_phase[0] for solution_phase in solution_phases], state['temperature'])
                
            stable_phases.append(solution_phases + pure_condensed_phases)            

        return stable_phases
    
    def _get_mole_fraction_element_by_phase(self):
        mole_frac_element_by_phase = { element: [] for element in self.elements }
        for element in self.elements:
            for state_index, state in enumerate(self.output.values()):
                mole_frac_element_by_phase[element].append([]) # Add an empty list corresponding to the given state
                for stable_phase, _ in self.stable_phases[state_index]: # [0] is the state name, [1] is the mole fraction
                    # First get phase type
                    if stable_phase in list(state['solution phases'].keys()):
                        phase_type = 'solution phases'
                    else:
                        phase_type = 'pure condensed phases'

                    # Then save mole fraction of element by phase for the given element and stable phase to the output dictionary
                    if element in list(state[phase_type][stable_phase]['elements'].keys()):
                        mole_frac_element_by_phase[element][state_index]\
                            .append( \
                                ( stable_phase, state[phase_type][stable_phase]['elements'][element]['mole fraction of element by phase']) \
                                    )
                    else:
                        mole_frac_element_by_phase[element][state_index].append( ( stable_phase, 0 ) )
        return mole_frac_element_by_phase



# ---------------------------------
#       Auxillary Functions
# ---------------------------------


def component_fractions_to_element_fractions(component_fractions, unique_elements):
    """Function for converting from a component-wise composition definition, where a compound
    is specified by mole fractions of components (e.g. NaCl-MgCl2-PuCl3 60-20-20 mol%) to an element
    wise definition (e.g. Na-Cl-Mg-Pu 0.6-1.6-0.2-0.2 mol)
    
    Parameters
    ----------
        unique_elements: A list of the string names of the unique components, e.g. 
                         unique_elements = ['Na', 'Cl', 'Mg', 'Pu']
        component_fractions: A dictionary of component string-mole fraction pairs, e.g. 
                             component_fractions = {'Na Cl': 0.6, 'Mg Cl_2': 0.2, 'Pu Cl_3': 0.2}
    
    Returns
    -------
        A list containing the mole fractions of each of the unique elements in the same order as unique_elements
    """
    element_fractions = [0]*len(unique_elements)

    for key, mol_frac in component_fractions.items():
        elements = [element.split('_')[0] for element in key.split()]
    
        # Now, get stochiometric coefficients of each element
        coefs = np.zeros(len(key.split()))
        for index, key in enumerate(key.split()):
            if len(key.split('_')) == 1:
                coefs[index] = 1
            else:
                coefs[index] = int(key.split('_')[1])
        
        # Now multiply mole fraction by stoichiometric coef to get total mole fraction of a given element in a component
        mol_frac_elements = mol_frac * coefs
        for i, element in enumerate(elements):
            element_index = unique_elements.index(element)
            element_fractions[element_index] += mol_frac_elements[i]

    # Note that, since the calculation only depends on the mole fractions of the elements, rather than absolute amounts
    # we may normalize the element fractions to 1
    return element_fractions


def get_unique_elements(components: list) -> list:
    """Gets the unique elements in a list of strings representing components"""
    elements = set()
    for component in components:
        elements.update( [ element.split('_')[0] for element in component.split() ] )

    return list(elements)



def solubility_calculation(temp: float, press: float, unit_ratio_of_other_components: dict, component_to_vary: str, n_comp_step: int, \
                           thermochimica_path: Path, output_path: Path, output_name: str, data_file: Path, \
                           compstart: float=0.0, compstop: float=1.0, script_name: str="thermoInput.ti", \
                           fuzzy: bool=False) -> thermoOut:
    """Function for performing a sequence of thermochimica calculations at a fixed temperature and pressure
    corresponding to varying one component of a system (e.g. PuCl3) while keeping the others in a fixed ratio
    
    Parameters:
    -----------
        temp: Fixed temperature in units of (K) to run the series of thermochimica calculations at
        press: Fixed pressure in units of (atm) to run the series of thermochimica calculations at
        unit_ratio_of_other_components: A dictionary containing the mole fractions of a unit solution of only
                                        the other components. The sum of the mole fractions MUST equal 1, because this
                                        serves only to define the ratio of these components to each other (which remains
                                        fixed, while the total mole fraction is varied as component_to_vary is varied)
                                        e.g. {'Na Cl': 0.8, 'Mg Cl_2': 0.2} corresponding to a ratio of 0.8:0.2 = 4:1
        component_to_vary: A string representing the component to vary e.g. 'Pu Cl_3'
        n_comp_setp: Number of composition steps to partition the composition interval [0,1] into
        thermochimica_path: Absolute path to the thermochimica directory. e.g. Path("/home/user/thermochimica")
        output_path: Path of the output directory (relative to the thermochimica/outputs directory). e.g. if the
                     thermochimica_path is given as above, Path("../../myProject/outputs") would write thermochimica outputs to
                     the 'outputs' subdirectory of the 'myProject' directory in the /home/user directory
        output_name: The name of the output file, e.g. \"output.json\"
        data_file: The path to the datafile (relative to the where this script is run) that is used to run the thermochimica calculation
                   e.g. Path("data/MSTDB-TC_V3.0_Chlorides_No_Functions_8-2.dat")
        compstart: The initial mole fraction of `component_to_be_varied`
        compend: The final mole fraction of `component_to_be_varied`
        script_name: Name of the input script to be created in the directory where this script is run from, by default \"runThermochimica\"
                     this will be deleted after each calculation
        fuzzy: Whether the calculation should be run with fuzzy stochiometry
    
    Returns:sum(element_fractions)
    --------
        A thermoOut object containing the results of the calculations in order as mol fraction of component_to_vary groes from 0 to 1
    """

    # Get unique elements and a list of components
    components = list(unit_ratio_of_other_components.keys()) + [component_to_vary]
    unique_elements = get_unique_elements(components)

    mol_fracs_component_to_vary = np.linspace(compstart, compstop, n_comp_step)

    # Perform the calculation for varying component_to_vary mole fractions
    output = thermoOut()
    for mol_frac_component_to_vary in mol_fracs_component_to_vary:
        # First get the mole fraction of each of the components
        frac_other_components = { component: unit_frac*(1-mol_frac_component_to_vary)  \
                                 for component, unit_frac in unit_ratio_of_other_components.items() }
        components = {**frac_other_components, component_to_vary: mol_frac_component_to_vary}

        # Write input script
        tstart = temp
        tend   = temp
        ntstep = 1
        pstart = press
        pend   = press
        npstep = 1
        masses = component_fractions_to_element_fractions(components, unique_elements)
        
        thermoTools.WriteInputScript(script_name, str(data_file), unique_elements, \
                                     tstart, tend, ntstep, pstart, pend, npstep, masses, fuzzyStoichiometry=fuzzy)

        # Run script
        thermoTools.RunInputScript(script_name, jsonName=str(output_path / output_name), thermochimica_path=str(thermochimica_path),
                                   noOutput=True)

        output_json_path = thermochimica_path / 'outputs' / output_path / output_name
        output.add_output(output_json_path)

    # Now remove run script and output
    os.remove(script_name)
    os.remove(str(output_path / output_name))

    return output