from lamost_cannon import *
l = full_sample()
#l = l.iloc[700000:].reset_index(drop=True)
normalize_spectra_loop(l)
#normalize_spectra_loop_parallel(l, 16)
