# MOlten Static Salt (MOSS) Reactor Senior Design Project
This repository contains all of the code necessary for reproducing the calculations done for the authors' senior design project for NRE 4351 at Georgia Tech.

## Table of Contents

- [Fuel/coolant salt optimization](./saltOptimization/README.md)
- [Core height scoping calculation](./coreHeightScoping/README.md)

## Requirements

Most of the scripts in this repository depend on having access to specific databases/scripts. The requirments are as follows:
- **The Molten Salt Thermal Properties Database (MSTDB)**: A database for calculating thermophysical/thermochemical properties of molten salts which is maintained by Oak Ridge National Laboratory (instructions for getting access are [here](https://mstdb.ornl.gov/about/))
  - All scripts that use the MSTDB for thermophsyical properties calculations reference a [symbolic link](#symlinks) in the root of the repository, which will have to be replaced with a valid symlink to _your_ clone of the MSTDB-TP repo.
- **Thermochimica**: A equilibrium thermodynamics solver which is used throughout for solubility calculations and worst-case corrosion calculations. The repository is hosted [here](https://github.com/ORNL-CEES/thermochimica).
  - The 

### Setup

#### Symlinks