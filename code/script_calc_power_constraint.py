from __future__ import division
from pylab import *
from scipy import interpolate

import blahut, helpers

#==============================================================================
#%% Load and preformat the data
#Derive 2-D Density with Linear Interpolation of Gaussian KDE
#==============================================================================
data = load('./../data/IEDM2014_PCM_Partial_Reset.npz')
helpers.dict2global(data)
Pyx, x, y = blahut.Q(V, R, nx=2000, ny=2000)

#==============================================================================
# %% Calculate Constrained Capacity
#==============================================================================
f = load('./../data/transistorFD12.npz')
Vg = f['Vg']
I = f['Isd'].T[-1]
I_spline = interpolate.UnivariateSpline(Vg, I)

Ix = I_spline(x)
Vsd = 4
t_reset = 300e-9 #seconds
energy = (Ix * Vsd * t_reset * 1e12) / 2 #pJ

#s go from 1e-2 to 1
n = 10
s_list = logspace(-4, -1, n)

func = lambda s: blahut.blahut_arimoto(Pyx=Pyx,
                                       tolerance=1e-3,
                                       iterations=10,
                                       e=energy,
                                       s=s,
                                       debug=True)
iters = [('s', s_list)]
outs = ['C','Px', 'E']
res = helpers.loop(func, iters, outs, savefile='./../npz/Cenergy_loop.npz')


