import pandas as pd
import numpy as np

def check_photometry(data):
    fltr = (data['G']!=data['G'])|(data['eG']!=data['eG'])
    if np.count_nonzero(fltr):
        data.loc[fltr, 'mag_use'] = \
            data.loc[fltr].applymap(lambda x: np.array(['J', 'H', 'K']))
    fltr2 = np.copy(fltr.values)
    fltr2 &= (data['K']!=data['K'])|(data['eK']!=data['eK'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['J', 'H']))
    fltr2 = np.copy(fltr)
    fltr2 &= (data['H']!=data['H'])|(data['eH']!=data['eH'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['J', 'K']))
    fltr2 = np.copy(fltr)
    fltr2 &= (data['J']!=data['J'])|(data['eJ']!=data['eJ'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['H', 'K']))
    fltr2 = ~np.copy(fltr)
    fltr2 &= (data['K']!=data['K'])|(data['eK']!=data['eK'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['J', 'H', 'G']))
    fltr2 = ~np.copy(fltr)
    fltr2 &= (data['H']!=data['H'])|(data['eH']!=data['eH'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['J', 'K', 'G']))
    fltr2 = ~np.copy(fltr)
    fltr2 &= (data['J']!=data['J'])|(data['eJ']!=data['eJ'])
    if np.count_nonzero(fltr2):
        data.loc[fltr2, 'mag_use'] = \
            data.loc[fltr2].applymap(lambda x: np.array(['H', 'K', 'G']))
    return data
