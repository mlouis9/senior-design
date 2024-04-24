# Average Coolant Case
This is the coolant calculation for the _average_ fuel tube, i.e. the power density is 100 kW/L, although we might expect the hottest fuel tube to have a power density nearly twice this.

## Processing Scripts
There are two paraview postprocessing scripts, [radialAverages.py](radialAverages.py), [innerSurfaceAverages.py](innerSurfaceAverages.py) which read the output (in this case from [613](613)) and write the extracted radially averaged, and clad surface averaged fuelds to csv files in the extracts directory.

## Run Script
There is a PBS script used for submitting this job on a high performance computing cluster. This script was written specifically to be run on Sawtooth (a cluster available via the [NCRC](https://inl.gov/ncrc/)).

## Supplementary Calculations
The supplementary calcualtions needed to determine the heat flux (used as a boundary condition), and postprocess the axial clad temperature profile (used as a boundary condition for the fuel calculation) is contained in [supplementaryCalcs.ipynb](supplementaryCalcs.ipynb).