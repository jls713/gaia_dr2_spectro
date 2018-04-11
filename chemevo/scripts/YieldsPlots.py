# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd

def yields_plot(fname,outname,masses=[1.,2.,5.],Z=0.001,lim=[10**-18,50.],label='',clabel=r'$M_\star/\mathrm{M}_\odot$'):

    data = pd.read_csv(fname,sep=r'\s+')
    data = data[data.Z==Z].reset_index(drop=True)

    fig = plt.figure(figsize=[15.,4.])
    norm = colors.Normalize(vmin=np.min(masses),vmax=np.max(masses))
    c_m = plt.cm.viridis
    s_m = plt.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    for m in masses:
        plt.plot(data[data.M==m].Yield.values,'.-',markersize=10,color=s_m.to_rgba(m))
    cmm=plt.colorbar(s_m)
    cmm.set_label(clabel)

    ele = data.Element.values
    indexes = np.unique(ele, return_index=True)[1]
    el= [ele[index] for index in sorted(indexes)]

    plt.xticks(range(len(el)),el)
    plt.semilogy()
    plt.ylabel(r'$M_\mathrm{el}/\mathrm{M}_\odot$')
    plt.ylim(lim[0],lim[1])
    plt.annotate(label,
                 xy=(0.95,0.05),horizontalalignment='right',verticalalignment='bottom',
                 xycoords='axes fraction',clip_on=False,fontsize=24)
    plt.tight_layout()
    plt.savefig(outname)

if __name__ == '__main__':
    yields_plot('tmp','output.png')
