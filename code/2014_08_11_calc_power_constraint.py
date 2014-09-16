from __future__ import division
from pylab import *
from scipy import interpolate

from myCodeCopy import JessePlot, blahut, CapacityTools, helpers
dict2global, loop = helpers.dict2global, helpers.loop
Q = CapacityTools.Q

#==============================================================================
#%% Load and preformat the data
#Derive 2-D Density with Linear Interpolation of Gaussian KDE
#==============================================================================
data = load('./../data/FD12W_2_PS_4.npz')
dict2global(data)

Pyx, x, y = Q(V, R, nx=2000, ny=2000)


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


#x = X[:, 0]
#y = Y[0, :]
#Pyx = np.rot90(Z)

#s go from 1e-2 to 1

n = 100

s_list = logspace(-4, -1, n)

C_list = [ 0 for i in range(n)]
Px_list =  [ 0 for i in range(n)]
E_list = [ 0 for i in range(n)]

for i in arange(n):
    s = s_list[i]

    C, Px, E = blahut.blahut_arimoto(Pyx,
                                     tolerance=1e-3,
                                     iterations=1000,
                                     e = energy,
                                     s = s,
                                     debug=True)

    C_list[i], Px_list[i], E_list[i] = C, Px, E




#%% Save Px
savez('./../npz/Px_power.npz', Px=Px)


#%% Load Px
Px = load('./../npz/Px_power.npz')['Px']


#%% Save Full Data
data = dict(C_list=C_list, Px_list=Px_list, E_list=E_list)
savez('./../npz/CEnergy.npz', data=data)

