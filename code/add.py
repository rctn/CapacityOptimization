#%%


import operator
import numpy as np
import helpers




def add_to_database(result, index, goal, savefile='test.npz'):
    ''' 
    Checks if result is better than database at index, and if so, replaces
    the database at index for each of outs and saves it to disk. Requires that all 
    dictionary values in savefile be of the same length (index is same).

    Parameters
    ----------
    result
        (dict)
    index
        (TUPLE)
    goal
        (tuple) ('variable','max' or 'min')
    outs
        (list of strings)
    savefile
        (string to .npz file)
    '''
    try:
        database = np.load(savefile)
    except IOError:
        print 'Wrong filename'
        return
        
    goal_var, maxmin = goal
    goal_cmp = [operator.lt, operator.gt][['min', 'max'].index(maxmin)]
    
    data_val = database[goal_var][index]
    res_val = result[goal_var]

    # Check if better
    if goal_cmp(res_val, data_val):

        # Copy the Dict (Necessary for new save for some reason)
        temp_dict = {}
        for k, v in database.items():
            temp_dict[k] = v

        # Overwrite data
        for out in result.keys():
            if out not in ['iters_dict', 'order']:
                temp_dict[out][index] = result[out]

        # Save
        np.savez(savefile, **temp_dict)



#%%

def test(x, y, z):
    a = x + y + z
    b = x*y*z
    return a, b
    
func = lambda x, y: test(x, y, z=1)
iters = [('x', range(3)), ('y', range(5))]
outs = ['a', 'b']
 
res = helpers.loop(func, iters, outs, savefile='test.npz')
helpers.dict2global(res)

#%%
r = {'a':11, 'b':12}
goal = ('a','max')
index = (0,0)
outs = ['a', 'b']

dataz = add_to_database(r, index, goal, 'test.npz')


#%% 
data = np.load('test.npz')
data.items()