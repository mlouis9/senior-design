from pyprojroot import here
import os
import re

""" This script is for setting up environment variables for use with pyprojroot. I.e. each script in this repository can have
from modules.paths import ROOT, THERMOCHIMICA_FLUORIDE_DATA, etc. instead of manually including relative paths which may change
if files are moved.

Author: Matthew Louis
Email: matthewlouis31@gmail.com
"""

def get_database_version(file_name):
    version_pattern = r"\d+\.\d+" # Pattern for finding version number of database
    match = re.search(version_pattern, file_name)
    if match:
        return float(match.group())
    else:
        return 0.0

# Now search through all files in the thermochimica data directory
def get_database_files(THERMOCHIMICA_DATA, chloride_database_string, fluoride_database_string):
    # Get the list of files in the thermochimica data directory
    files = os.listdir(str(THERMOCHIMICA_DATA))

    # Iterate over the files and find the ones with the database strings
    chloride_matching_files = []
    fluoride_matching_files = []

    for file_name in files:
        if chloride_database_string in file_name:
            version = get_database_version(file_name)
            chloride_matching_files.append((version, file_name))
        elif fluoride_database_string in file_name:
            version = get_database_version(file_name)
            fluoride_matching_files.append((version, file_name))
        
    # If no files were found, raise an errror
    if not chloride_matching_files:
        raise ValueError(f"There are no data files containing the string {chloride_database_string} in the thermochimica data directory."
                        "Either the chloride database string is specified incorrectly in the paths.py script, or you have not "
                        "copied the MSTDB data files to the thermochimica data directory.")

    if not fluoride_matching_files:
        raise ValueError(f"There are no data files containing the string {fluoride_matching_files} in the thermochimica data directory."
                        "Either the fluoride database string is specified incorrectly in the paths.py script, or you have not "
                        "copied the MSTDB data files to the thermochimica data directory.")

    # Now choose the most recent file for each database
    chloride_database = max(chloride_matching_files, key=lambda x: x[0])[1]
    fluoride_database = max(fluoride_matching_files, key=lambda x: x[0])[1]
    return chloride_database, fluoride_database

# ----------------------------------------------------------
# Setting up environment variables for use with pyprojroot
# ----------------------------------------------------------

ROOT = here()
THERMOCHIMICA = here("thermochimica")
MSTDB_TP = here("mstdb-tp")
THERMOCHIMICA_DATA = here("thermochimica/data")

# Now find the files corresponding to the most recent versions of MSTDB-TC for the fluorides and chlorides

chloride_database_string = "Chlorides_No_Functions"
fluoride_database_string = "Fluorides_No_Functions"

chloride_database, fluoride_database = get_database_files(THERMOCHIMICA_DATA, chloride_database_string, fluoride_database_string)
THERMOCHIMICA_CHLORIDE_DATA = here(f"thermochimica/data/{chloride_database}")
THERMOCHIMICA_FLUORIDE_DATA = here(f"thermochimica/data/{fluoride_database}")

MSTDB_TP_DATA = here("mstdb-tp/Molten_Salt_Thermophysical_Properties.csv")
MSTDB_TP_RK_DATA = here("mstdb-tp/Molten_Salt_Thermophysical_Properties_RK.csv")
SALT_OPTIMIZATION = here("saltOptimization")

# ----------
# CFD Cases
# ----------

AVERAGE_COOLANT_CASE = here("cfd/3Dpincell/averageCoolant")
SHORT_FUEL_CASE = here("cfd/3Dpincell/shortFuel")
AVERAGE_FUEL_CASE = here("cfd/3Dpincell/averageFuel")