import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
import os
sys.path.append('.')
from dyneq import models, greybox

def plot2D(yref, yident, pl_results, figsize=(20,5), pl_option='vpl', n_ppl=4):
    fig = plt.figure(figsize=figsize)
    ax = plt.axes()

    z_max = np.empty(len(pl_results))
    z_min = np.empty(len(pl_results))
    t_z = np.empty(len(pl_results))

    time = yref.index.values
    for idx, pl in enumerate(pl_results):
        if pl_option == 'vpl':
            indices = (pl.vpl < 3.8)
            z_max[idx] = pl.z[indices].max()
            z_min[idx] = pl.z[indices].min()
        elif pl_option == 'ppl':
            indices = (pl.ppl < 3.8)
            z_max[idx] = pl.zhat[indices].max()
            z_min[idx] = pl.zhat[indices].min()
        t_z[idx] = yref.index.values[pl.idx]

    ax.fill_between(t_z, z_max, z_min, alpha=0.2, facecolor='darkseagreen')
    ax.plot(yident['tf.y'],'darkseagreen')
    ax.plot(yref['tf.y'],'darkred')
    plt.xlim([0,10])
    plt.xlabel('Time (s)')
    plt.ylabel('Active Power (kW)')
    return fig, ax

with open('results.p','rb') as file:
    res = pickle.load(file)

with open('yref.p','rb') as file:
    yref = pickle.load(file)

with open('yident.p','rb') as file:
    yident = pickle.load(file)
    

plot2D(yref, yident, res, figsize = (15,5),pl_option='ppl')
plt.show()
print('ha')