from __future__ import division
from pylab import *
from scipy import optimize

from myCodeCopy import JessePlot, blahut, CapacityTools, helpers
dict2global = helpers.dict2global
Q = CapacityTools.Q

#==============================================================================
#%% Load and preformat the data
#==============================================================================
data = load('./../data/FD12W_2_PS_4.npz')
dict2global(data)
Pyx, x, y = Q(V, R, nx=2000, ny=2000)

#==============================================================================
# %% Discrete inputs
#==============================================================================
reload(blahut)

new_run = False # Erases old file Only 1st Time
new_x0 = True

totsize = 100
nx = 30
x_range = (min(x),max(x))


C_in_temp = [0 for i in range(totsize)]
x_in_temp = [0 for i in range(totsize)]
Px_in_temp = [0 for i in range(totsize)]


for repeat in arange(100):
    for ix in arange(5, 13):
    #    if new_x0 or new_run:
    #        xinputs = rand(ix) * abs(x_range[1]-x_range[0]) + x_range[0]
    #        xinputs.sort()
    #        x0 = xinputs
    #    else:
    #        x0 = x_in[ix]
        x0 = np.random.choice(x0_max, ix, replace=False)
    
            
        bounds = [x_range for i in arange(ix)]
        params = (Pyx, x, y, ix)
        
         
        minimizer_kwargs = dict(method='L-BFGS-B',
                                args=params,
                                bounds=bounds)
            
        
        def objective(val, Pyx=Pyx, x=x, y=y, ix=ix):
            xinputs = val
            Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs, ydividers=None)
            C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                          tolerance=1e-4, 
                                          iterations = 100)
                                    
            if C > C_in_temp[ix]:        
                C_in_temp[ix] = C
                x_in_temp[ix] = xinputs
                Px_in_temp[ix] = Px
                data = dict(C_in_temp=C_in_temp, x_in_temp=x_in_temp, Px_in_temp=Px_in_temp)
                savez('./../npz/basinhopping_xin_temp.npz', data=data)    
    
    
            print '\nC:', C
            print 'xin:', xinputs
            print 'nx:', ix
            return -C
        
        
        
    
        result = optimize.basinhopping(objective, 
                                       x0=x0, 
                                       minimizer_kwargs=minimizer_kwargs,
                                       niter=500,
                                       stepsize=0.5,
                                       T=0.5)
    
        
        print result


#==============================================================================
# %% Combine temp XIN
#==============================================================================

r = load('./../npz/basinhopping_xin_temp.npz')['data'].item()
C_opt, x_in_opt, Px_opt = (r['C_in_temp'], 
                     r['x_in_temp'], 
                     r['Px_in_temp'])

r = load('./../npz/basinhopping.npz')['data'].item()
C_out, y_div, Px_out, C_in, x_in, Px_in = (r['C_out'], 
                                           r['y_div'], 
                                           r['Px_out'],
                                           r['C_in'], 
                                           r['x_in'], 
                                           r['Px_in'])
                     
#%%                     
                     
for i in arange(len(C_out)):                   
    if C_opt[i] > C_in[i]:
        C_in[i] = C_opt[i]
        x_in[i] = x_in_opt[i]
        Px_in[i] = Px_opt[i]

#%%

data = dict(C_out=C_out, C_in=C_in, x_in=x_in, y_div=y_div, Px_in=Px_in, Px_out=Px_out)
savez('./../npz/basinhopping.npz', data=data)    



#==============================================================================
# %% Discrete outputs
#==============================================================================
#
reload(blahut)
new_run = True # Erases old file Only 1st Time
new_x0 = True

totsize = 100
ny = 30
y_range = (min(y),max(y))

C_out_temp = [0 for i in range(totsize)]
y_div_temp = [0 for i in range(totsize)]
Px_out_temp = [0 for i in range(totsize)]

for niteration in arange(1000):
    for iy in [31, 63]:#arange(2, ny+1):
        if new_x0 or new_run:
            ydividers = rand(iy) * abs(y_range[1]-y_range[0]) + y_range[0]
            ydividers.sort()
            x0 = ydividers
        else:
            x0 = y_div[iy]
    
            
        bounds = [y_range for i in arange(iy)]
        params = (Pyx, x, y, iy)
        
         
        minimizer_kwargs = dict(method='L-BFGS-B',
                                args=params,
                                bounds=bounds)
    
        
        def objective(val, Pyx=Pyx, x=x, y=y, iy=iy):
            global C_out_temp, y_div_temp, Px_out_temp
    
            ydividers = array(val)
            ydividers.sort()
            Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs=None, ydividers=ydividers)
            C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                          tolerance=1e-4, 
                                          iterations = 100)
                                    
            print '\nC:', C
            print 'ydividers:', ydividers
            print 'ny:', iy
            
            if C > C_out_temp[iy]:        
                C_out_temp[iy] = C
                y_div_temp[iy] = ydividers
                Px_out_temp[iy] = Px
                data = dict(C_out_temp=C_out_temp, y_div_temp=y_div_temp, Px_out_temp=Px_out_temp)
                savez('./../npz/basinhopping_yout_temp.npz', data=data)    
            
            return -C
        
        
        def accept_test(f_new, x_new, f_old, x_old):
            '''
            Make it so there aren't two dividers with the same value
            '''
            if len(x_new) == len(set(x_new)):
                return True
            else:
                return False
        
    
        result = optimize.basinhopping(objective, 
                                       x0=x0, 
                                       minimizer_kwargs=minimizer_kwargs,
                                       niter=500,
                                       stepsize=2.0,
                                       T=0.5,
                                       accept_test=accept_test)
        
        
        
        print result


   

#==============================================================================
# %% Combine temp YOUT
#==============================================================================

r = load('./../npz/basinhopping_yout_temp.npz')['data'].item()
C_opt, y_div_opt, Px_opt = (r['C_out_temp'], 
                     r['y_div_temp'], 
                     r['Px_out_temp'])

r = load('./../npz/basinhopping.npz')['data'].item()
C_out, y_div, Px_out, C_in, x_in, Px_in = (r['C_out'], 
                                           r['y_div'], 
                                           r['Px_out'],
                                           r['C_in'], 
                                           r['x_in'], 
                                           r['Px_in'])
                
#%%
close('all'); plot(C_opt); plot(C_out); show()               
                
#%%                     
                     
for i in arange(len(C_out)):                   
    print C_out[i], C_opt[i]
    if C_opt[i] > C_out[i]:
        C_out[i] = C_opt[i]
        y_div[i] = y_div_opt[i]
        Px_out[i] = Px_opt[i]

#%%

data = dict(C_out=C_out, C_in=C_in, x_in=x_in, y_div=y_div, Px_in=Px_in, Px_out=Px_out)
savez('./../npz/basinhopping.npz', data=data)    





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




   


