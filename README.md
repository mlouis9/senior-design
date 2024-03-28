# MOlten Static Salt (MOSS) Reactor Senior Design Project
This repository contains all of the code necessary for reproducing the calculations done for the authors' senior design project for NRE 4351 at Georgia Tech.

## Table of Contents

- [coreHeightScoping](coreHeightScoping/README.md)
- [saltOptimization](saltOptimization/README.md)
  - [benchmarking](saltOptimization/benchmarking/README.md)
  - [corrosionCalcs](saltOptimization/corrosionCalcs/README.md)
  - [thermophysicalProperties](saltOptimization/thermophysicalProperties/README.md)
    - [benchmarking](saltOptimization/thermophysicalProperties/benchmarking/README.md)

## Requirements

Most of the scripts in this repository depend on having access to specific databases/scripts. The requirments are as follows:
- **The Molten Salt Thermal Properties Database (MSTDB)**: A database for calculating thermophysical/thermochemical properties of molten salts which is maintained by Oak Ridge National Laboratory (instructions for getting access are [here](https://mstdb.ornl.gov/about/))
  - All scripts that use the MSTDB for thermophsyical properties calculations reference a [symbolic link](#updating-the-symlinks) in the root of the repository, which will have to be replaced with a valid symlink to _your_ clone of the MSTDB-TP repo.
- **Thermochimica**: A equilibrium thermodynamics solver which is used throughout for solubility calculations and worst-case corrosion calculations. The official repository is hosted [here](https://github.com/ORNL-CEES/thermochimica). A fork of the repository containing updated python functions which are _necessary_ for the code in this repository to run is [here](git@github.com:mlouis9/thermochimica.git). When cloning this repository, run
  ```
  git checkout extra-options 
  ```
  as this repository is not updated with the official repository, it may be better to clone the official thermochimica repo, and then copy the `python` directory from my personal fork.

## Setup
After obtaining access to the [required codes/data](#requirements), some setup is still required before being able to run the code in this repository.
### Required Python Packages
First, install all of the required python packages, which are listed in the [requirements.txt](./requirements.txt) file. If pip is your python package manager, this can be done by running
```
pip install -r requirements.txt
```
in the root of the repository. If you're using mamba or conda as your package manager, run
```
<package-manager> install --file requirements.txt
```

### Updating the Symlinks
For the relative paths in this repository to work, symbolic linlks (aka symlinks, more information [here](https://en.wikipedia.org/wiki/Symbolic_link)) must be created in the root of this repository for both the `thermochimica` and `mstdb-tp` directories. The current symlinks must be deleted first, which can be done (in Linux) by running
```
rm thermochimica mstdb-tp
```
To create new symlinks (in Linux) run the following commands in the root of this repoisitory
```
ln -s <path-to-thermochimica-directory> ./thermochimica
ln -s <path-to-mstdb-tp-directory> ./mstdb-tp
```
NOTE: The thermochimica directory should be named `thermochimica`, and the thermophysical properties directory should be named `mstdb-tp`.

### Setting up Thermochimica
To run thermochimica calculations using the MSTDB thermochimical database, the relevant chem sage `.dat` files must be copied to the thermochimica data directory. To do this, you can run the following command (in Linux):
```
cp <path-to-mstdb>/Models\ and\ Documentation/*.dat <path-to-thermochimica-directory>/data
```

### Adding Additional Modules to PYTHONPATH
Finally, it is necessary to add the python scripts in the [modules](./modules/) directory, as well as the python utility script for running thermochimica calculations to your PYTHONPATH, otherwise the various import statements will fail. To do this, you may add the following lines to your `.bashrc` (in Linux)
```
export PYTHONPATH=<path-to-thermochimica>/python:$PYTHONPATH
export PYTHONPATH=<path-to-this-repo>/modules:$PYTHONPATH
```