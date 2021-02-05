import numpy as np
from scipy.io import loadmat
import os
from pathlib import Path
# from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


# plotting parameters
sns.set(font_scale=1.1)
sns.set_context("talk")
sns.set_palette(['#701f57', '#ad1759', '#e13342', '#f37651'])
transparent = False
markers = ['o','^','s']
plot_types = ['box','point']
estimator = np.median
est ='median'

# Number of microphones fixed, vary number of loudspeakers

path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(Path(__file__).parent,'figures')
if not os.path.exists(outpath):
    os.makedirs(outpath)

Ks = range(6,12) # number of loudspeakers
M = 12 # number of microphones
runs = 200

thresh = 10**-3 # lower bound for error

algs = ['SDR + LM', 'Wang']
suffs = ['real_data_clean_sdr_lm', 'wang_real_data_clean']

res_data = pd.DataFrame() # will load matlab results into DataFrame
all_err = np.zeros((runs,len(Ks),len(suffs)))
for mi,K in enumerate(Ks):

    N = M + K      # number of points
    print(mi,M,K)
      
    for suffi,suff in enumerate(suffs):
        filename_mat = os.path.join(suff,'matlab_%s_M%s_K%s.mat'%(suff,M,K))
        mat = loadmat(os.path.join(path,filename_mat))
        err_aft = np.real(mat['err_lm'])
        
        all_err[:,mi,suffi] = err_aft[:,0]
        print('minimum',algs[suffi],M,K, np.min(all_err[:,mi,suffi]))

        res = pd.DataFrame(err_aft[:,0])
        res['M'] = M
        res['K'] = K
        res.rename(columns={0: 'err'}, inplace=True)
        res = res.assign(Algorithm=algs[suffi])
        
        res_data = pd.concat([res_data,res],axis=0)

# res_data.clip(lower=thresh)
for plot_type in plot_types:
    args ={}
    if plot_type=='point':
        args['markers'] = markers
        args['s'] = 0.75

    ax = sns.catplot(data=res_data, x='K', y='err', hue='Algorithm', linewidth=1, dodge=True, kind=plot_type, legend_out=True, legend=False, aspect=1.3,estimator=estimator,**args)
    ax.set_ylabels("Localization Error [m]")
    plt.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.1))
    plt.tight_layout()
    plt.savefig(os.path.join(outpath,"loc_error_compare_%s_%s_%splot.png"%(est,suff,plot_type)), transparent=transparent)
