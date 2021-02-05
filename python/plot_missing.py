import numpy as np
from scipy.io import loadmat
import os
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
sns.set(font_scale=1.1)
sns.set_context("talk")
sns.set_palette(['#701f57', '#ad1759', '#e13342', '#f37651'])


transparent = False
num_col = 3 # for plotting
linestyles = ['-']
markers = ['o','^','s']
plot_types = ['box','point']
estimator = np.median
est ='median'

path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(Path(__file__).parent,'figures')
if not os.path.exists(outpath):
    os.makedirs(outpath)


noises = [1e-3,1e-4,1e-5,1e-6,0] # noise stdev
noiselabels =  [r'$10^{-3}$',r'$10^{-4}$',r'$10^{-5}$',r'$10^{-6}$',r'$0$']

deletions = np.linspace(0.01,0.3,10)

thresh = 10**-3 # lower bound for error

# list of figures to generate
figs = {}
algs = ['Complete Data','Missing Data']
suffs = ['bef_lm', 'missing']
Ms = [12] # number of microphones

figs[0] = [algs,suffs,Ms]

algs = ['Missing Data']
suffs = ['missing']
Ms = [8,10,12] # number of microphones
figs[1] = [algs,suffs,Ms]

for key in figs.keys():
    algs = figs[key][0]
    suffs = figs[key][1]
    Ms = figs[key][2]
    
    res_data = pd.DataFrame()
    for di,dele in enumerate(deletions):
        if di<1 or di>6:
            continue
        print(di,dele)
        
        for M in Ms:

            K = M          # number of sources
            N = M + K      # number of points

            for suffi,suff in enumerate(suffs):
                filename_mat = os.path.join(suff,'matlab_%s_M%s.mat'%(suff,M))
                mat = loadmat(os.path.join(path,filename_mat))
                err_bef = np.real(mat['err_bef']) # loc error before refinement
                err_aft = np.real(mat['err_lm']) # loc error after refinement

                if suff=='missing':
                    err_aft = err_aft[...,di]
                
                res = pd.DataFrame(err_aft, columns=noises).reset_index().melt(id_vars = 'index')\
                .rename(columns={'index':'runs'})
                res.rename(columns={"variable": "sigma", "value": "err"}, inplace=True)
                res = res.assign(Algorithm=algs[suffi])
                res = res.assign(M=M)
                res = res.assign(Missing='%.2f%%'%(dele*100))

                res_data = pd.concat([res_data,res])

        
    
    for plot_type in plot_types:
        args ={}
        if plot_type=='point':
            args['markers'] = markers
            args['s'] = 0.75
            args['linewidth'] = 1
        
        if len(Ms)>1:  
            ax = sns.catplot(data=res_data, x='sigma',y='err', col='Missing', hue='M', order=noises,dodge=True,kind=plot_type,legend_out=True,col_wrap=num_col,legend=False,estimator=estimator,**args)
        else:
            ax = sns.catplot(data=res_data, x='sigma',y='err', col='Missing', hue='Algorithm', order=noises,dodge=True,kind=plot_type,legend_out=True,col_wrap=num_col,legend=False,estimator=estimator,**args)
        
        ax.set(xticklabels=noiselabels)
        ax.set(yscale='log')
        ax.set_xlabels("Noise $\sigma$")
        ax.set_ylabels("Localization Error [m]")
        ax.set(ylim=[thresh,10**1])
        ax.axes[0].legend(loc='lower left')
        plt.tight_layout()
        plt.savefig(os.path.join(outpath,"loc_error_compare_%s_%s_%s_%splot.png"%(est,suff,key,plot_type)), transparent=transparent)
