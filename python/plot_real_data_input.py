import numpy as np
import os, datetime
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy.io import loadmat
import pandas as pd
sns.set(font_scale=1.1)
sns.set_palette('rocket')

path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(Path(__file__).parent,'figures')
if not os.path.exists(outpath):
    os.makedirs(outpath)


M = 12
K = 65

d = 3       # ambient dimension
c = 343        # speed of sound

filename_mat = 'real_data_M%s_K%s.mat'%(M,K)
mat = loadmat(os.path.join(path,filename_mat))
T = mat['T']
W = mat['W']
X_r = mat['X_r']

fig, axes = plt.subplots(1, 1)
out0 = axes.matshow(T, aspect='auto')
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(out0, cax=cax) 
axes.grid(False)    

fig.tight_layout()

plt.savefig(os.path.join(outpath,"real_data_T.png"))

fig, axes = plt.subplots(1, 1)
out0 = axes.matshow(W, aspect='auto')
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(out0, cax=cax) 
axes.grid(False)    

fig.tight_layout()

plt.savefig(os.path.join(outpath,"real_data_W.png"))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X_r[0,:], X_r[1,:], X_r[2,:], marker='o',color='b')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
M = X_r.shape[1]
for mi in range(M):
    ax.text(X_r[0,mi], X_r[1,mi], X_r[2,mi], mi+1)
fig.tight_layout()

plt.savefig(os.path.join(outpath,"real_data_Xr.png"))


