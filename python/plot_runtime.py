import numpy as np
from scipy.io import loadmat
import os
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

# plotting parameters
sns.set(font_scale=1.1)
sns.set_context("talk")
sns.set_palette(['#701f57', '#ad1759', '#e13342', '#f37651'])
transparent = False
markers = ['o','^','s']

# Plot runtime as we increase the number of sensors
path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(Path(__file__).parent,'figures')
if not os.path.exists(outpath):
    os.makedirs(outpath)


Ms = range(7,13)

res_data = pd.DataFrame()

suff = 'runtime_comp_err_noono'
algs = ['SDR + LM', 'Wang']
filename_mat = os.path.join(suff,'matlab_%s_M12.mat'%(suff))
mat = loadmat(os.path.join(path,filename_mat))

ours = np.real(mat['runtime_ours']).T 
wang = np.real(mat['runtime_wang']).T

print(ours.shape)

res = pd.DataFrame(ours,columns=Ms).reset_index().melt(id_vars = 'index')\
            .rename(columns={'index':'runs'})
res.rename(columns={"variable": "M", "value": "runtime"}, inplace=True)
res = res.assign(Algorithm=algs[0])

res_data = pd.concat([res_data,res],axis=0)

res = pd.DataFrame(wang,columns=Ms).reset_index().melt(id_vars = 'index')\
            .rename(columns={'index':'runs'})
res.rename(columns={"variable": "M", "value": "runtime"}, inplace=True)
res = res.assign(Algorithm=algs[1])

res_data = pd.concat([res_data,res],axis=0)

print(res_data)

ax = sns.catplot(data=res_data, x='M',y='runtime',hue='Algorithm', linewidth=1,dodge=True, kind='point',aspect=1.3,legend=False,scale=0.75,markers=markers)#,legend_out=True)
ax.set_xlabels("M=K")
ax.set_ylabels("seconds")
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(outpath,"runtime_compare_%s_pointplot.png"%suff), transparent=transparent)
