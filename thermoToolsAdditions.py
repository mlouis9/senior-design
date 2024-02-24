import thermoTools
import json
from pathlib import Path
import numpy as np
import os
import copy
import pseudoBinaryPhaseDiagramFunctions as pbpd
import matplotlib.pyplot as plt
import math
from scipy.spatial import ConvexHull

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
            self.temperatures = dict()
            bad_keys = []
            for key, state in self.output.items():
                try:
                    self.temperatures.update({key: state['temperature']})
                except:
                    # There is a problem with this state, i.e. it is either completely empty or missing required data, print
                    # the key for debugging purposes
                    print(f"Warning: state {key} is missing data! Excluding this state from output")
                    
                    # Now addd this to the list of bad keys
                    bad_keys.append(key)

            # Now remove all bad keys
            for key in bad_keys:
                self.output.pop(key)

            self.elements = list( set( element for state in  list(self.output.values()) \
                                 for element in list(state['elements'].keys()) ) ) # Get unique elements corresponding to each of the
                                                                                   # states, should be the SAME for all states (although
                                                                                   # different mol fractions are possible). NOTE that
                                                                                   # thermochimica does not add elements with zero mole
                                                                                   # fractions to output.json, though they should still
                                                                                   # be included for our purposes
            self.stable_phases = self._get_stable_phases()
            self.mole_fraction_element_by_phase = self._get_mole_fraction_element_by_phase()

        else:
            self.temperatures = dict()
            self.stable_phases = dict()
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


    def _get_stable_phases(self):
        # Below this tolerance, set phase fraction = 0
        phase_include_tol = 1e-8

        stable_phases = dict()
        for state_name, state in self.output.items():
            phase_types = ['solution phases', 'pure condensed phases']
            phases = []
            for phase_type in phase_types:
                for name, phase in state[phase_type].items():
                    if phase['moles'] > phase_include_tol:
                        # First add phase composition
                        phase_composition = []
                        for element in self.elements:
                            # NOTE not all elements are in all phases
                            if element in phase['elements']:
                                phase_composition.append( phase['elements'][element]['mole fraction of phase by element'] )
                        phases.append( (name, phase['moles'], phase_composition) )
            stable_phases.update({state_name: phases})

        return stable_phases
    
    def _get_mole_fraction_element_by_phase(self):
        mole_frac_element_by_phase = { element: {} for element in self.elements }
        for element in self.elements:
            for state_index, state in self.output.items():
                mole_frac_element_by_phase[element].update({state_index: []}) # Add an empty list corresponding to the given state
                for stable_phase, *_ in self.stable_phases[state_index]: # [0] is the phase name, [1] is the mole fraction, [2] is the composition
                    # First get phase type
                    if stable_phase in list(state['solution phases'].keys()):
                        phase_type = 'solution phases'
                    else:
                        phase_type = 'pure condensed phases'

                    # Then save mole fraction of element by phase for the given
                    # element and stable phase to the output dictionary
                    if element in list(state[phase_type][stable_phase]['elements'].keys()):
                        mole_frac_element_by_phase[element][state_index]\
                            .append( \
                                ( stable_phase, state[phase_type][stable_phase]['elements'][element]['mole fraction of element by phase']) \
                                    )
                    else:
                        mole_frac_element_by_phase[element][state_index].append( ( stable_phase, 0 ) )
        return mole_frac_element_by_phase


class pseudoBinaryDiagram(thermoOut):
    """A class which extends thermoOut for use in making pseudo binary phase diagrams"""
    def __init__(self, left_endmember_composition: dict, right_endmember_composition: dict, out_file: Path = None):
        # First create a thermoOut object with the output file
        super().__init__(out_file)
        self.left_endmember = left_endmember_composition
        self.right_endmember = right_endmember_composition

        # Get mass labels for plotting
        left_endmember_masses = component_fractions_to_element_fractions(left_endmember_composition, self.elements)
        right_endmember_masses = component_fractions_to_element_fractions(right_endmember_composition, self.elements)
        self.mass_labels = get_mass_labels(left_endmember_masses, right_endmember_masses, self.elements)

        # Now add a new attribute that is a dictionary of the mole fraction of right endmember and state
        self.mol_frac_right_endmember = {\
            state_key: \
                element_to_component_fractions_pseudo_binary(\
                    left_endmember_composition, right_endmember_composition, \
                    { element_name: element_dict['moles'] for element_name, element_dict in self.output[state_key]['elements'].items() }\
                                                                )\
                                          for state_key, _ in self.output.items() }
        
        # Now get regions
        self._get_phase_regions()

    # def _get_boundary_points(self):
    #     # Assuming a common carrier element (e.g. Cl or F), the number of components that determine the number of degrees of freedom is given
    #     # by the number of elements -1 (the carrier), thus the phase boundaries occur when F = 1 = 2 + C - P ==> P = 1 + C = number of elements

    #     self.phase_boundaries = [ phases for phases in self.stable_phases if len(phases) == len(self.elements)]
    def _get_phase_regions(self):
        """Function for getting the unique phase regions for each set of stable phases, these will partition the composition/temperature space
        into regions, after which we can find the boundaries"""

        # Assuming a common carrier element (e.g. Cl or F), the number of components that determine the number of components is given
        # by the number of elements -1 (the carrier), thus the phase boundaries occur when F = 1 = 2 + C - P ==> P = 1 + C = number of elements.
        # So the regions with two degrees of freedom (areas on the phase diagram) are those for which P = number of elements - 1

        region_points = { state_key: phases for state_key, phases in self.stable_phases.items() if len(phases) == len(self.elements) - 1 }
        self.regions = dict()

        for state_key, phases in region_points.items():
            phase_names = set(phase[0] for phase in phases)

            # As a dictionary key, we use a frozenset which is an immutable set object which will allow us to use an order-indepenent
            # way of accessing the points associated with a given set of phases
            phases_key = frozenset(phase_names)

            if phases_key in self.regions:
                self.regions[phases_key].append( ( self.mol_frac_right_endmember[state_key], self.temperatures[state_key] ) )
            else:
                self.regions.update({phases_key: [(self.mol_frac_right_endmember[state_key], self.temperatures[state_key])] })

        # Now convert to numpy arrays
        self.regions = {key: np.array(region) for key, region in self.regions.items()}


    def _filter_phase_points(self, threshold=0.1):
        """This function serves to eliminate points in phase regions that are far from the centroid of the rest of the points - which
        in general are the result of roundoff errors and numerical uncertainties, and can ruin our plots"""

        # We want to decrease the area of this convex hull (by excluding outliers) until the fractional decrease in area is less than a certain
        # threshold. Note this process filters out outliers because excluding them results in a large fractional area decrease of the convex hull

        def fractional_area_decrease(points, hull, point):
            """Calculate the fractional decrease in area of the convex hull by excluding a point."""
            original_area = hull.volume
            mask = np.any(points != point, axis=1)
            reduced_hull = ConvexHull(points[mask])
            reduced_area = reduced_hull.volume
            return (original_area - reduced_area) / original_area

        for region_key, points in self.regions.items():
            # Compute the convex hull of the points in a given phase region. NOTE the convex hull of a set of points is the smallest polygon that
            # contains all of the points

            if points.shape[0] < 4:
                # Cannot make a convex hull with less than 3 points, and cannot remove a point (and still calculate area) without 4 points
                continue

            hull = ConvexHull(points)
            filtered_points = points.copy()

            # We iterate over the boundary of the convex hull, and see which point (upon exclusion) results in the maximal decrease in the convex hull
            # area, then exclude this point and repeat until the fractional decrease is less than a threshold

            while True:
                boundary_vertices = filtered_points[hull.vertices]
                fractional_area_decreases= []
                for boundary_vertex in boundary_vertices:
                    fractional_area_decreases.append((boundary_vertex, fractional_area_decrease(filtered_points, hull, boundary_vertex)))

                # Find point that maximizes the fractional area decrease
                max_point, max_decrease = max(fractional_area_decreases, key=lambda x: x[1])
                if max_decrease < threshold:
                    break

                # Now filter out point
                mask = np.any(filtered_points != max_point, axis=1)
                filtered_points = filtered_points[mask]

                # Now create a new convex hull
                if filtered_points.shape[0] > 3:
                    hull = ConvexHull(filtered_points)
                else:
                    # Cannot calculate the area of a convex hull with two points
                    break
               
            # Now overwrite region points with filtered points
            self.regions[region_key] = filtered_points


    def plot_phase_regions(self, plot_marker='.'):
        # Initialize the phase_region_plot attribute
        self.phase_region_plot = plotObject()
        self.phase_region_plot.fig = None
        self.phase_region_plot.ax = None

        self.phase_region_plot.fig, self.phase_region_plot.ax = plt.subplots()

        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.regions.keys()))))
        for phase_names, phase_region in self.regions.items():
            c = next(color)

            x_points = phase_region[:, 0]
            y_points = phase_region[:, 1]

            self.phase_region_plot.ax.plot(x_points, y_points, plot_marker, c=c, label = '+'.join(phase_names))

        self.phase_region_plot.ax.set_xlim(0,1)
        title = " $-$ ".join(self.mass_labels)
        self.phase_region_plot.ax.set_title(r'{0} phase diagram'.format(title))
        self.phase_region_plot.ax.set_xlabel(r'Mole fraction {0}'.format(self.mass_labels[1]))
        self.phase_region_plot.ax.set_ylabel("Temperature [K]")
        self.phase_region_plot.ax.legend(loc='upper right', bbox_to_anchor=(1.5, 1))

        plt.show()


class plotObject:
    """An empty class for holding phase region plots"""
    def __init__(self) -> None:
        pass

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


def elements_from_component_key(component_key):
    return [element.split('_')[0] for element in component_key.split()]


def element_to_component_fractions_pseudo_binary(left_endmember_composition, right_endmember_composition, element_composition):
    """Function for taking element fractions (read from thermochimica output.json) and translating them to component fractions
    endmembers for a pseudo binary system. This may result in errors if the incorrect endmembers are provided, because there may
    be no solutions"""

    # First, get the element compositions of the endmembers and normalize them (thermochimica results are normalized for pesudo binary
    # calculaltions)
    unique_elements = []
    for key in left_endmember_composition.keys():
        unique_elements += elements_from_component_key(key)
    for key in right_endmember_composition.keys():
        unique_elements += elements_from_component_key(key)
    unique_elements = list(set(unique_elements))
    left_endmember_element_composition = component_fractions_to_element_fractions(left_endmember_composition, unique_elements)
    left_endmember_element_composition = list(left_endmember_element_composition / sum(left_endmember_element_composition))

    right_endmember_element_composition = component_fractions_to_element_fractions(right_endmember_composition, unique_elements)
    right_endmember_element_composition = list(right_endmember_element_composition / sum(right_endmember_element_composition))

    # Now, find the element(s) in the right endmember that are not present in the left-endmember
    for index, _ in enumerate(unique_elements):
        if (left_endmember_element_composition[index] == 0) and (right_endmember_element_composition[index] != 0):
            unique_endmember_index = index
            unique_endmember_element = unique_elements[index]
            break

    # Now can easily calculate mole fraction of the right endmember via
    if unique_endmember_element in element_composition:
        frac_right_endmember = element_composition[unique_endmember_element]/right_endmember_element_composition[unique_endmember_index]
    else:
        frac_right_endmember = 0 # Thermochimica doesn't print elements with zero mole fractions in output.json


    # Now, make sure that endmembers were properly specified. If not, the mole fraction computed above may not result in the correct element_composition
    calculated_element_composition = { element: 0.0 for element in unique_elements }

    left_endmember_elements = []
    right_endmember_elements = []
    for component_key in left_endmember_composition.keys():
        left_endmember_elements += elements_from_component_key(component_key)
    for component_key in right_endmember_composition.keys():
        right_endmember_elements += elements_from_component_key(component_key)

    for index, element in enumerate(unique_elements):
        if element in left_endmember_elements:
            calculated_element_composition[element] += (1-frac_right_endmember)*left_endmember_element_composition[index]
        if element in right_endmember_elements:
            calculated_element_composition[element] += frac_right_endmember*right_endmember_element_composition[index]
    
    # Now get rid of zero values (which won't be present in teh original 'element_composition')
    calculated_element_composition = {key: value for key, value in calculated_element_composition.items() if value != 0}

    def are_dicts_almost_equal(a, b, tolerance=1e-9):
        # Check if the keys are the same
        if set(a.keys()) != set(b.keys()):
            return False
        
        # Check if the values are almost equal for each key
        for key in a:
            if not math.isclose(a[key], b[key], abs_tol=tolerance):
                return False
        
        return True

    # Now convert to element fractions to compare with the original 'element_composition'
    # Note, there will be some roundoff error when computing the "calculated_element_composition" dict, so we assert that these dictionaries must
    # be within some absolute tolerance (by default taken to be 1e-9) of each other
    assert are_dicts_almost_equal(calculated_element_composition, element_composition), f"Inconsistent endmembers {left_endmember_composition} and {right_endmember_composition} specified for the given element_composition {element_composition}! No solution can be found."

    return frac_right_endmember


def get_mass_labels(left_endmember_masses: dict, right_endmember_masses: dict, elements_used: list) -> list:
    """Utility for getting the labels of endmembers in a pseudo binary phase diagram from their compositions. This
    code is more or less ripped directly from pseudoBinaryPhaseDiagramGUI.py
    
    Parameters:
    -----------
        left_endmember_masses: A dictionary containing component identifier - mole fraction pairs representing the left
                               endmember
        right_endmember_masses: A dictionary containing component identifier - mole fraction pairs representing the right
                               endmember
        elements_used: A list containing the string identifiers of the unique elements used in both of the endmembers

    Returns:
    --------
        A list of two strings, containing the proper pylatex commands for rendering the endmember molecular formulas (with
        subscripts!)
    """
    mass_labels = ['', '']
    n_elements_used = len(elements_used)

    # If halogens are present, make sure they are the last element to appear in each label
    # by convention. This can be done by moving them to the end of the elements_used list
    # and also permuting the right and left endmember masses lists correspondingly
    if ('Cl' in elements_used) ^ ('F' in elements_used): # Note we use xor `^` to exclude albeit unphsical cases where BOTH F and Cl are present
        if 'Cl' in elements_used:
            halogen_index = elements_used.index('Cl')
        else:
            halogen_index = elements_used.index('F')

        # First copy the input list to avoid side effects
        lists_to_be_permuted = copy.deepcopy([elements_used, left_endmember_masses, right_endmember_masses])
        for list_to_be_permuted in lists_to_be_permuted:
            # Then shuffle the elements
            halogen = list_to_be_permuted.pop(halogen_index)
            list_to_be_permuted.append(halogen)

        # Now save modified lists under the proper names
        elements_used = lists_to_be_permuted[0]
        left_endmember_masses = lists_to_be_permuted[1]
        right_endmember_masses = lists_to_be_permuted[2]

    for i in range(n_elements_used):
        if left_endmember_masses[i] > 0:
            mass_labels[0] += elements_used[i]
            if left_endmember_masses[i] != 1:
                if int(left_endmember_masses[i]) == left_endmember_masses[i]:
                    mass_labels[0] += f'$_{ {int(left_endmember_masses[i])} }$'
                else:
                    mass_labels[0] += f'$_{ {left_endmember_masses[i]} }$'
        if right_endmember_masses[i] > 0:
            mass_labels[1] += elements_used[i]
            if right_endmember_masses[i] != 1:
                if int(right_endmember_masses[i]) == right_endmember_masses[i]:
                    mass_labels[1] += f'$_{ {int(right_endmember_masses[i])} }$'
                else:
                    mass_labels[1] += f'$_{ {right_endmember_masses[i]} }$'
    return mass_labels


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

def pseudo_binary_calculation(thermochimica_path: Path, output_path: Path, output_name: str, data_file: Path, xlo: float, xhi: float, nxstep: int, \
                              tlo: float, thi: float, ntstep: int, elements_used: list, left_endmember_composition: dict, \
                              right_endmember_composition: dict, press: float=1, tunit: str='K', punit: str='atm', munit: str='moles', \
                              input_file_name: str='runThermochimica.ti', fuzzy: bool=True) -> pbpd.diagram:
    """Function for performing a pseudo binary calculation. This performs most of the tedious IO tasks, and properly converts from component compositions
    to element fractions, etc.
    
    Parameters:
    -----------
        thermochimica_path:
        output_path:
        output_name:
        data_file:
        xlo: The low composition (usually 0.0)
        xhi: The high composition (usually 1.0)
        nxstep: The number of composition steps to partition the interval [xlo, xhi] into
        tlo: The low temperature of the phase diagram
        thi: The high temperature of the phase diagram
        ntstep: The number of temperature steps to partition the interval [tlo, thi] into
        elements_used: A list containing the string identifiers of the unique elments which the endmembers are composed of
        left_endmember_composition: A dictionary containing the component identifiers and mole fractions of the left endmember
        right_endmember_composition: A dictionary containing the component identifiers and mole fractions of the right endmember
        press:
        tunit:
        punit:
        munit:
        input_file_name:
        fuzzy:
    
    Returns:
    --------
        A pseudo binary phase diagram object containing the relevant calculation results

    """
    left_endmember_masses = component_fractions_to_element_fractions(left_endmember_composition, elements_used)
    right_endmember_masses = component_fractions_to_element_fractions(right_endmember_composition, elements_used)

    sum1 = sum(left_endmember_masses)
    sum2 = sum(right_endmember_masses)

    mass_labels = get_mass_labels(left_endmember_masses, right_endmember_masses, elements_used)

    # Now normalize left and right endmember masses, NOTE this is REQUIRED for the plotting and postprocessing to behave
    left_endmember_masses = list( left_endmember_masses / sum1 )
    right_endmember_masses = list( right_endmember_masses / sum2 )

    plane = [left_endmember_masses, right_endmember_masses]

    if tunit == 'K':
        tshift = 0
    elif tunit == 'C':
        tshift = 273.15
    else:
        assert(f"Invalid temperature unit \'{tunit}\' provided")

    mint = tlo + tshift
    maxt = thi + tshift

    calc = pbpd.diagram(data_file, active=True, interactivePlot=True, inputFileName="runThermochimica.ti", \
                    outputFileName=str(output_path / output_name), thermochimicaPath=thermochimica_path)

    calc.initRun(press, tunit, punit, plane, sum1, sum2, mint, maxt, elements_used, mass_labels, munit, tshift, fuzzy=fuzzy)
    calc.runCalc(xlo, xhi, nxstep, tlo, thi, ntstep)
    calc.processPhaseDiagramData()
    return calc