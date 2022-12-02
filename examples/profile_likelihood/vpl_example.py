'''
Example of VPL estimation using a PVD1.
'''

if __name__ == '__main__':
    import os
    import sys
    sys.path.append('.') # Project root
    sys.path.append('./examples')
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    #plt.style.use('science')

    from dyneq import models, greybox, vpl
    import model_builder

    ########################## GETTING THE MODELS ##########################
    refmodel, equivmodel = model_builder.build(from_executables=False)

    ###################### DEFINING THE IDENTIFICATOR ######################

    # Variables to compute residuals
    ref_vars = ['source.P', 'source.Q']
    equiv_vars = ['source.P', 'source.Q']
    
    # Simulating the reference measurements
    yref = refmodel.simulate()

    # Define the parameters and map them to the positions in the list
    param_map = {
        0: ('solar.controller.Ddn'),
        1: ('solar.controller.dqdv'),
        2: ('solar.controller.fdbd'),
        3: ('solar.controller.Xc'),
        4: ('solar.controller.fr_recov'),
        5: ('solar.controller.v0'),
        6: ('solar.controller.v1'),
        7: ('solar.controller.Vt0'),
        8: ('solar.controller.Vt1'),
        9: ('solar.controller.Ft0'),
        10:('solar.controller.Ft1'),
        11:('solar.controller.Ft2'),
        12:('solar.controller.Ft3'),
        13: ('solar.controller.Tg'), 
    }

    # Define lower and upper bounds
    lb = np.array([ 0,   0,   -0.1,  0, 0, 0.9,  1,  0, 0, 0, 0, 1, 1, 0.015])
    ub = np.array([0.1, 0.1, 0.0001, 1, 1,  1,  1.1, 1, 1, 1, 1, 2, 2, 0.05])

    identificator = greybox.LeastSquares(model=equivmodel, measurements=yref,
                                        ref_vars=ref_vars, target_vars=equiv_vars,
                                        param_map=param_map)

    # Setting the weights
    #weights = 0.01*np.ones(len(yref))
    #idxs = np.logical_and(yref.index > 3.28,yref.index < 8)
    #weights[idxs] = 1
    #identificator.weights = np.append(weights, weights)

    # Assuming theta_hat is the optimum
    theta_hat = [0.04,1e-3, -0.001, 0, 1, 0.95, 1.05, 0.9, 0.91,
                 0.998, 0.999, 1.005, 1.01, 0.02]

    ##################### START OF PL ESTIMATION ###########################
    ident_opts = {'method':'trf', 'tr_solver':'exact','verbose':2, 'bounds':(lb,ub), 'jac':'2-point','max_nfev':150}

    var_names = ['source.P']*1
    indices = [345]

    # CHOOSE HERE BETWEEN MULTICORE AND SINGLECORE
    opt = vpl.PLOptions(theta_hat, 1, 0.05, multicore=False, 
                        save_intermediate=True, ident_opts=ident_opts)
    
    
    vplest = vpl.VPLEstimator(equivmodel, identificator)
    
    start = time.perf_counter()
    res = vplest.multiple_estimates(var_names, indices, opt)
    end = time.perf_counter()
    print(f'Estimation time: {end-start} seconds')
    
    plt.plot(yref['source.P'].values, 'x')
    for idx in indices:
        plt.plot(idx, res[idx].z.min(),'rx')
        plt.plot(idx, res[idx].z.max(),'gx')
    plt.savefig('dummy.png')