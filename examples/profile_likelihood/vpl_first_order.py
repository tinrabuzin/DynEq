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
    import pickle

    from dyneq import models, greybox,pl

    # Defining paths
    build_path =  os.path.abspath('./data/tmp').replace("\\","/")
    model_path = os.path.abspath('./data/models/Modelica/DynEq/package.mo').replace("\\","/")

    output_specs=['tf.y']

    experiment_settings = {
        'startTime':0, 'stopTime':10, 'numberOfIntervals':int(3/0.1), 'method':'"rungekutta"'
    }
    
    
    
    # Creating the reference model
    model = models.OMModel(existing = False, model_file_path = model_path, model_name = 'DynEq.Equivalents.first_order', 
                                binary_folder = build_path.replace("\\","/"), output_specs=output_specs,
                                experiment_settings=experiment_settings)

    yref = model.simulate()
    
    # Variables to compute residuals
    ref_vars = ['tf.y']
    equiv_vars = ['tf.y']

    # Define the parameters and map them to the positions in the list
    param_map = {
        #0: ('tf.k'), # Current controller time constant
        0: ('tf.T'),
    }

    
    identificator = greybox.LeastSquares(model=model, measurements=yref,
                                        ref_vars=ref_vars, target_vars=equiv_vars,
                                        param_map=param_map)
    identificator.weights = np.ones(yref['tf.y'].values.shape)/0.1
    #lb = np.array([0,0.01])
    #ub = np.array([100,10])
    lb = np.array([0.01])
    ub = np.array([10])
    ident_opts = {'x0':[0.2], 'method':'trf', 'tr_solver':'exact', 'bounds':(lb,ub), 'jac':'3-point','gtol':1e-12}
    res = identificator.identify(**ident_opts)
    # Selection of trajectory and indices on that trajectory
    indices = [*range(0,len(yref), 1)]
    var_names = ['tf.y']*len(indices)
    
    opt = pl.PLOptions(res.x, 0.1, 0.05, multicore=True, 
                        save_intermediate=False, ident_opts=ident_opts, 
                        stop_criteria='vpl', stepping_method='past_based')

    plest = pl.PLEstimator(model, identificator)
    plres = plest.multiple_estimates(var_names, indices, opt)

    with open('results.p', 'wb') as file:
        pickle.dump(plres, file)
