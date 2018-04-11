# gaia_dr2_prep

Distances, ages and masses for spectroscopic surveys with Gaia

Datasets

1. APOGEE
2. LAMOST
3. GALAH
4. RAVE
5. GES

All the paths are specified in config.json

1. I download isochrones to the dir/isochrones folder in config.json using get_isochrones.py
2. I construct a dense extinction map with reddening.py (that is saved in dir/extinction folder in config.json)

The code (for action finding) needs my tact code.
