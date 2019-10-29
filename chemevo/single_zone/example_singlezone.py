import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../scripts/')
from chemevo import *

f = chem_evo_data('output.hdf5')

f.plot_time('Fe')
plt.show(block=True)
