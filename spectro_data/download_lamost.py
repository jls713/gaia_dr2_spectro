from lamost_cannon import *
l = full_sample()
#l=l.iloc[].reset_index(drop=True)
download_spectra_loop(l)
#download_spectra_loop_parallel(l, 16)
