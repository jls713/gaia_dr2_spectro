from lamost_cannon import *
l = full_sample().reset_index(drop=True)
#l=l.iloc[679389:].reset_index(drop=True)
#download_spectra_loop(l)
download_spectra_loop_parallel(l, 200)
