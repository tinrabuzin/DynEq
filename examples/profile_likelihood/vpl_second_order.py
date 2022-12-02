'''
Example of VPL estimation using a simple first order model.
'''
  
if __name__ == '__main__':
    import os
    import sys
    import copy
    sys.path.append('.') # Project root
    import numpy as np
    import matplotlib.pyplot as plt
    import time

    from dyneq import models, greybox,vpl

    # Defining paths
    build_path =  os.path.abspath('./data/tmp').replace("\\","/")
    model_path = os.path.abspath('./data/models/Modelica/DynEq/package.mo').replace("\\","/")

    output_specs=['tf.y']

    experiment_settings = {
        'startTime':0, 'stopTime':10, 'numberOfIntervals':int(10/0.1), 'method':'"rungekutta"'
    }
    
    
    
    # Creating the reference model
    model = models.OMModel(existing = False, model_file_path = model_path, model_name = 'DynEq.second_order', 
                                binary_folder = build_path.replace("\\","/"), output_specs=output_specs,
                                experiment_settings=experiment_settings)

    yref = model.simulate()
    
    # Variables to compute residuals
    ref_vars = ['tf.y']
    equiv_vars = ['tf.y']

    # Define the parameters and map them to the positions in the list
    param_map = {
        0: ('tf.k'), # Current controller time constant
        1: ('tf.w'),
        2: ('tf.D')
    }

    
    identificator = greybox.LeastSquares(model=model, measurements=yref,
                                        ref_vars=ref_vars, target_vars=equiv_vars,
                                        param_map=param_map)
    identificator.weights = np.ones(yref['tf.y'].values.shape)/0.1
    lb = np.array([0,0.1, 0.05])
    ub = np.array([100,10, 2])
    ident_opts = {'x0':[1,3,0.5], 'method':'trf', 'tr_solver':'exact', 'bounds':(lb,ub), 'jac':'3-point','gtol':1e-12}

    # Selection of trajectory and indices on that trajectory
    indices = [*range(0,101)]
    var_names = ['tf.y']*len(indices)
    
    
    # CHOOSE HERE BETWEEN MULTICORE AND SINGLECORE
    opt = vpl.VPLOptions([1, 3, 0.5], 0.01, 0.05, multicore=True, 
                        save_intermediate=True, ident_opts=ident_opts)

    vplest = vpl.VPLEstimator(model, identificator)
    start = time.perf_counter()
    res = vplest.multiple_estimates(var_names, indices, opt)
    end = time.perf_counter()
    print(f'Estimation time: {end-start} seconds')
    plt.plot(yref['tf.y'].values, 'x')
    for idx in indices:
        plt.plot(idx, res[idx].z.min(),'rx')
        plt.plot(idx, res[idx].z.max(),'gx')
    plt.savefig('dummy.png')
