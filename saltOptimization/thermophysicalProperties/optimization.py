from pathlib import Path
from thermophysicalProperties import Database
import matplotlib.pyplot as plt
import numpy as np
from frozendict import frozendict
import pint
ureg = pint.UnitRegistry(auto_reduce_dimensions=True)
from uncertainties import ufloat
import os
from scipy.misc import derivative

script_dir = Path(os.path.abspath('')).resolve().parent

mstdb_tp_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties.csv')
mstdb_tp_rk_path = Path('/home/mlouis9/mstdb-tp/Molten_Salt_Thermophysical_Properties_2.1.0_RK.csv')

db = Database(mstdb_tp_path, mstdb_tp_rk_path)