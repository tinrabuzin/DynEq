import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
import os
sys.path.append('.')
sys.path.append('./examples')
import pickle
import numpy as np
from dyneq import vpl
import matplotlib.tri as mtri
#plt.style.use('science')

def plot2D(yref, yident, pl_results, figsize=(12,8), pl='ppl', n_ppl=4):
    fig = plt.figure()
    ax = plt.axes()

    z_max = np.empty(len(pl_results))
    z_min = np.empty(len(pl_results))
    t_z = np.empty(len(pl_results))

    time = yref.index.values
    for idx, pl in enumerate(pl_results):
        z_max[idx] = pl.z.max()
        z_min[idx] = pl.z.min()
        t_z[idx] = yref.index.values[pl.idx]

    ax.fill_between(t_z, z_max, z_min, alpha=0.2, facecolor='darkseagreen')
    ax.plot(yident['source.P'],'darkseagreen')
    ax.plot(yref['source.P'],'darkred')

    return fig, ax


def plot3D(pl_results, figsize=(12,8), pl='ppl',n_ppl=4):
    fig = plt.figure()
    ax = Axes3D(fig)
    # Make background transparent
    ax.patch.set_alpha(0.0)
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))


    z_max = np.empty(len(pl_results))
    z_min = np.empty(len(pl_results))
    indices = np.array([*range(325,805,5)])
    for idx, pl in enumerate(pl_results):
        z_max[idx] = pl.z.max()
        z_min[idx] = pl.z.min()

    # Plot the limits
    #ax.plot(indices, z_max,'darkseagreen')
    #ax.plot(indices, z_min, 'darkseagreen')
    
    # Shade between min and max
    x = np.append(indices,indices[::-1])
    y = np.append(z_max, z_min[::-1])
    z = np.zeros(x.shape)
    verts =  [[*zip(x,y,z)]]
    ax.add_collection3d(Poly3DCollection(verts,facecolor='darkseagreen', alpha=0.5))
    
    # Ploting profile likelihoods
    for pl in pl_results[::n_ppl]:
        x = np.ones(pl.z.shape)*pl.idx
        y = pl.z
        if pl == 'ppl':
            z = pl.ppl
        else:
            z = pl.chi2
        ax.plot3D(x, y, z, 'grey', alpha=0.5)
    ax.set_zlim([0,5])
    ax.set_ylim([14.5,18.5])
    ax.set_xlim(indices.min(),indices.max())
    return fig, ax


indices = [*range(325,805,5)]
var_names = ['source.P']*len(indices)
pickle_files = [f'vpl_DynEq.Equivalent_ZIP_source.P_{idx}.p' for idx, var_name in zip(indices, var_names)]


plres = []
for pickle_file in pickle_files:
    with open(pickle_file, 'rb') as file:
        plres.append(pickle.load(file))

with open('yref.p', 'rb') as file:
    yref = pickle.load(file)
with open('y.p', 'rb') as file:
    y = pickle.load(file)


plot2D(yref, y, plres, pl='ppl')
plt.show()