import pandas as pd
import sys
sys.path.append('../')
from cross_match import crossmatch_2mass_ids_gaia_version2
g = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input.hdf5')
g2 = crossmatch_2mass_ids_gaia_version2(g, dist_max=10.)
