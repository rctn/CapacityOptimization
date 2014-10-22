from __future import division
from pylab import *

import operator, itertools




    
def objective(val, Pyx=Pyx, x=x, y=y, ix=ix):
    xinputs = val
    Pyx_sub, x_sub, y_sub = blahut.quantize(Pyx, x, y, xinputs, ydividers=None)
    C, Px = blahut.blahut_arimoto(Pyx_sub, 
                                  tolerance=1e-4, 
                                  iterations=100)
    return -C
    
    

minimizer_kwargs = dict(method='L-BFGS-B',
                        args=params,
                        bounds=bounds)

search_kwargs = dict(x0=x0, 
                     minimizer_kwargs=minimizer_kwargs,
                     niter=500,
                     stepsize=0.5,
                     T=0.5)


optimizer = lambda obj: optimize.basinhopping(obj, **search_kwargs)


x0 =[0]
outs = ['C', 'x', 'Px']
goal = {'C':'max'}



def func():
    x0 = np.random.choice(x0_max, ix, replace=False)

        
    bounds = [x_range for i in arange(ix)]
    params = (Pyx, x, y, ix)
    




def search(x0, goal, iters, outs, objective_func, savefile='test.npz', search_kwargs={}, n_trials=1):

    if not search_kwargs:
        search_kwargs = dict(minimizer_kwargs=minimizer_kwargs,
                             niter=500,
                             stepsize=0.5,
                             T=0.5)
    search_kwargs['x0']=x0



    goal_var, maxmin = goal.items()
    goal_cmp = [operator.lt, operator.gt][['min', 'max'].index(maxmin)]

    
    #Make sure iters is a list of tuples
    if len(iters[0]) < 2:
        iters = [iters]
        
    iter_keys, iter_values = zip(*iters)
    n_iters = [len(v) for v in iter_values]
#    values = itertools.product(*iter_values)

    #initialize n-dim list
    empty_list = ''
    for n in n_iters[::-1]:
        empty_list = [empty_list] * n


    #Load Database
    try:
        database = load(savefile)
    except IOError:
        database = {}
        for out in outs:
            if out:
                database[out] = empty_list            




    def add_to_database(result, index,
                        database=database, goal_var=goal_var, 
                        goal_cmp=goal_cmp, savefile=savefile, outs=outs):
        ''' 
        Checks if result is better than database at index, and if so, replaces
        the database and saves it to disk.
        '''
        data_val = database[goal_var][index]
        res_val = result[goal_var]

        if goal_cmp(res_val, data_val):
            
            for out in outs:                
                database[out][index] = result[out]
                
            print '%.2e better than %.2e at index %d' % (result, data_val, index)
            np.savez(savefile, **database)
            


    for repeat in arange(n_trials): 
                
        
        
        def new_objective(x, index=index, goal_var=goal_var):
            result = objective_func(x)
            add_to_database(result, index)
            return result
        

        result = optimize.basinhopping(new_objective, **search_kwargs)
        