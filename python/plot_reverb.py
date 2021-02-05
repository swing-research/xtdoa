import numpy as np
from scipy.io import loadmat
import os
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
sns.set(font_scale=1.1)
sns.set_palette(['#701f57', '#ad1759', '#e13342', '#f37651'])

transparent = False
linestyles = ['-']
markers = ['o','^','s','X']
plot_types = ['box','point']
estimator = np.median
est ='median'

path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(Path(__file__).parent,'figures')
if not os.path.exists(outpath):
    os.makedirs(outpath)

M = 12 # number of microphones
thresh = 10**-3 # lower bound for error

figs = {}

algs = ['Two-stage','SDR+LM']
suffs = ['reverb_aux_wang','reverb_aux']

res_data = pd.DataFrame()       
for suffi,suff in enumerate(suffs):
    filename_mat = os.path.join(suff,'matlab_%s_M%s.mat'%(suff,M))
    mat = loadmat(os.path.join(path,filename_mat))
    err_aft = np.real(mat['err_lm']) # loc error
    err_tdoa = np.real(mat['err_tdoa']) # tdoa error
    rt60 = mat['est_reverb'] # estimated rt60
    rt_rounded = np.round(rt60,decimals=1)
    print(rt_rounded)
    
    res = pd.DataFrame({'err':err_aft.reshape(-1,),'RT':rt60.reshape(-1,),'err_tdoa':err_tdoa.reshape(-1,),'RT60':rt_rounded.reshape(-1,)})
    res = res.assign(Algorithm=algs[suffi])
    res = res.assign(M=M)
    print(res.head())
    
    res_data = pd.concat([res_data,res])


res_data['err'] = res_data['err'].clip(lower=thresh)
aspect = 1.5
for plot_type in plot_types:
    args ={}
    if plot_type=='point':
        args['markers'] = markers
        args['s'] = 0.75
    ax = sns.catplot(data=res_data, x='RT60', y='err', col='M', hue='Algorithm', linewidth=1, dodge=True, kind=plot_type, legend_out=True, col_wrap=1, legend=False,estimator=estimator,aspect=aspect,**args)
    ax.set(yscale='log')
    ax.set_xlabels("Reverberation time [s]")
    ax.set_ylabels("Localization Error [m]")
    ax.set(ylim=[thresh,10**1])
    ax.axes[0].legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(outpath,"loc_error_reverb.png"), transparent=transparent)

    ax = sns.catplot(data=res_data[res_data.Algorithm.eq(algs[-1])], x='RT60', y='err_tdoa', col='M', hue='Algorithm', linewidth=1, dodge=True, kind=plot_type, legend_out=True, col_wrap=1, legend=False,estimator=estimator,aspect=aspect,**args)
    ax.set(yscale='log')
    ax.set_xlabels("Reverberation time [s]")
    ax.set_ylabels("TDOA Estimation Error [s]")
    plt.tight_layout()
    plt.savefig(os.path.join(outpath,"tdoa_error_reverb.png"), transparent=transparent)