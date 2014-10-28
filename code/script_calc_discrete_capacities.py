from __future__ import division
from pylab import *
from scipy import optimize

import JessePlot, blahut, helpers

#==============================================================================
#%% Load and preformat the data
#==============================================================================
data = load('./../data/IEDM2014_PCM_Partial_Reset.npz')
V = data['V']
R = log10(data['R'])
Pyx, x, y = blahut.Q(V, R, nx=2000, ny=2000)

#==============================================================================
# %% Discrete inputs
#==============================================================================

filename = './../npz/IEDM/optimal_discrete_inputs.npz'

x_range = (min(x),max(x))
nbits = arange(1, 9)

#Initialize the File
try:
    data = load(filename)
except IOError:
    totsize = nbits.size
    for var in ['C_in', 'C_out']:
        locals()[var] = zeros(totsize)
    for var in ['x_in', 'y_div', 'Px_in', 'Px_out']:
        locals()[var] = [[0 for j in range(i)] for i in range(totsize)]
    data = dict(C_in=C_in, x_in=x_in, Px_in=Px_in, nbits=nbits)
    savez(filename, **data)    



new_x0 = True


#Search
for repeat in arange(100):
    idx = np.random.choice(nbits) - 1
    states = 2**(idx + 1)

    if new_x0 or new_run:
        x0 = rand(states) * abs(x_range[1]-x_range[0]) + x_range[0]
        x0.sort()
    else:
        x0 = data['x_in'][idx]
    
        
    bounds = [x_range for i in arange(states)]
    params = (Pyx, x, y, idx)
    
     
    minimizer_kwargs = dict(method='L-BFGS-B',
                            args=params,
                            bounds=bounds)
        
    
    def objective(val, Pyx=Pyx, x=x, y=y, idx=idx):
        xinputs = val
        xinputs.sort()
        
        Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs, ydividers=None)
        C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                      tolerance=1e-4, 
                                      iterations = 100)
                                
        result = {'C_in':C, 'Px_in':Px, 'x_in':xinputs}
        helpers.add_to_database( result=result, index=idx, 
                                 goal=('C_in', 'max'), filename=filename )

        print '\nC:', C
        print 'x_in:', xinputs
        print 'idx:', idx

        return -C
    

    result = optimize.basinhopping(objective, 
                                   x0=x0, 
                                   minimizer_kwargs=minimizer_kwargs,
                                   niter=100,
                                   stepsize=0.5,
                                   T=0.5)
    
    print result







#==============================================================================
# %% Equally spaced inputs / outputs Optimized
#==============================================================================

reload(blahut)

nbits=8
C_equal = zeros((nbits, nbits))
C_equal_bounds = zeros((nbits, nbits, 4)) #xmin, xmax, ymin, ymax
states = 2**(arange(nbits)+1)
xy = meshgrid(states, states)

xr = [min(x), (min(x) + max(x))/2, max(x)]
yr = [min(y), (min(y) + max(y))/2, max(y)]


def objective(val, Pyx=Pyx, x=x, y=y, nx=nx, ny=ny):
    global C_equal, C_equal_bounds
    
    xmin, xmax, ymin, ymax = val   

    xinputs = linspace(xmin, xmax, nx)                 
    ydividers = linspace(ymin, ymax, ny+1)
    ydividers = ydividers[1:-1]  
        
    Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs, ydividers)
    C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                  tolerance=1e-7, 
                                  iterations = 100)

    
    if C > C_equal[ix, iy]:        
        C_equal[ix, iy] = C 
        C_equal_bounds[ix, iy, :] = (xmin, xmax, ymin, ymax) 
        data = dict(C_equal=C_equal,
                    C_equal_bounds=C_equal_bounds)
        savez('./../npz/C_equal_optimized.npz', data=data)    
                                    

    
    print '\nC:', C
    print 'nx:', nx, xmin, xmax
    print 'ny:', ny, ymin, ymax

    return -C


def accept_test(f_new, x_new, f_old, x_old):
    '''
    Make it so there aren't two dividers with the same value
    '''
    if len(x_new) == len(set(x_new)):
        return True
    else:
        return False



#%%

for repeats in arange(1000):
    for ix in arange(nbits):
        for iy in arange(nbits):
            ny, nx = (xy[0][ix,iy], xy[1][ix,iy])
    
            params = (Pyx, x, y, nx, ny)
            bounds = [ (xr[0], xr[1]), (xr[1], xr[2]), (yr[0], yr[1]), (yr[1], yr[2]) ]
    
            minimizer_kwargs = dict(method='L-BFGS-B',
                                    args=params,
                                    bounds=bounds)
    
            x0 = array([xr[0], xr[2], yr[0], yr[2]])
        
            result = optimize.basinhopping(objective, 
                                           x0=x0, 
                                           minimizer_kwargs=minimizer_kwargs,
                                           niter=500,
                                           stepsize=1.0,
                                           T=0.5,
                                           accept_test=accept_test)
        
        
    
#%%
close('all')
matshow(C_equal)
show()
    

#%%    
data = dict(C_equal=C_equal,
            C_equal_bounds=C_equal_bounds)
savez('./../npz/C_equal_optimized.npz', data=data)    

#%%
r = load('./../npz/C_equal_optimized.npz')['data'].item()
C_equal= (r['C_equal'])
C_equal_bounds= (r['C_equal_bounds'])



#==============================================================================
# %% Nonuniform Spacing
#==============================================================================


reload(blahut)

nbits=8
C_nonequal = zeros((nbits, nbits))
x_in = [ [0 for i in range(nbits)] for j in range(nbits)]
y_out = [ [0 for i in range(nbits)] for j in range(nbits)]
Px_in = [ [0 for i in range(nbits)] for j in range(nbits)]


states = 2**(arange(nbits)+1)

xy = meshgrid(states, states)

xr = [min(x), (min(x) + max(x))/2, max(x)]
yr = [min(y), (min(y) + max(y))/2, max(y)]


def objective(val, Pyx=Pyx, x=x, y=y, nx=nx, ny=ny):
    global C_nonequal, x_in, y_out, Px_in
    
    xinputs = val[0:nx]
    ydividers = val[nx:]
    
    xinputs.sort()
    ydividers.sort()
        
    Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs, ydividers)
    C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                  tolerance=1e-7, 
                                  iterations = 1000)
    
    if C > C_nonequal[ix, iy]:        
        C_nonequal[ix, iy] = C 
        x_in[ix][iy] = xinputs
        y_out[ix][iy] = ydividers
        Px_in[ix][iy] = Px

        data = dict(C_nonequal=C_nonequal,
                    x_in = x_in, 
                    y_out = y_out,
                    Px_in = Px_in)
        savez('./../npz/C_nonequal_optimized.npz', data=data)    
                                    
    
    print '\nC:', C
    print 'x:', nx, xinputs
    print 'y:', ny, ydividers

    return -C


def accept_test(f_new, x_new, f_old, x_old):
    '''
    Make it so there aren't two dividers with the same value
    '''
    if len(x_new) == len(set(x_new)):
        return True
    else:
        return False



#%%

def find_peaks(Px):
    peak = logical_and( logical_and(Px > append(0, Px[:-1]), Px > append(Px[1:],0)), Px > 2e-4)
    return where(peak)


Px0 = load('./../npz/CEnergy.npz')['data'].item()['Px_list'][0]
x0_max = x[find_peaks(Px0)]

y0_max = load('./../npz/basinhopping.npz')['data'].item()['y_div'][10]




for it in range(10000):
    for ix, iy in [ (1,1), (2,2) ]:
        ny, nx = (xy[0][ix,iy], xy[1][ix,iy])
    
#        xinputs = rand(nx) * abs(xr[2]-xr[0]) + xr[0]
#        ydividers = rand(ny) * abs(yr[2]-yr[0]) + yr[0]    
#        xinputs = sort(xinputs)
#        ydividers = sort(ydividers)

#        xinputs = linspace(xmin, xmax, nx)                 
#        ydividers = linspace(ymin, ymax, ny+1)
#        ydividers = ydividers[1:-1]  

        xinputs = np.random.choice(x0_max, nx, replace=False)
        ydividers = np.random.choice(y0_max, ny-1, replace=False)
            

        x0 = np.append(xinputs, ydividers)
    
        params = (Pyx, x, y, nx, ny)
        bounds = [ (xr[0], xr[2]) for i in arange(nx) ]
        bounds.extend( [(yr[0], yr[2]) for i in arange(ny-1)] )
    
        minimizer_kwargs = dict(method='L-BFGS-B',
                                args=params,
                                bounds=bounds)
    
        result = optimize.basinhopping(objective, 
                                       x0=x0, 
                                       minimizer_kwargs=minimizer_kwargs,
                                       niter=50,
                                       stepsize=0.5,
                                       T=0.5,
                                       accept_test=accept_test)
            
            
    
#%%
close('all')
matshow(C_equal)
show()
    

#%%    
data = dict(C_nonequal=C_nonequal,
            x_in = x_in, 
            y_out = y_out,
            Px_in = Px_in)
savez('./../npz/C_nonequal_optimized.npz', data=data)  

#%%
r = load('./../npz//C_nonequal_optimized.npz')['data'].item()
for name in ['C_nonequal','x_in','y_out','Px_in']:
    globals()[name] = r[name]




   


