'''
Parameter identification for Modelica models with differential evoluton. 
'''
import sys
sys.path.append('.')
import numpy as np
import matplotlib.pyplot as plt
from dyneq import models, greybox, vpl
import model_builder

########################## GETTING THE MODELS ##########################
refmodel, equivmodel = model_builder.build(from_executables=True)

###################### DEFINING THE IDENTIFICATOR ######################

# Simulating the reference model
yref = refmodel.simulate()

# Variables to compute residuals
ref_vars = ['source.P', 'source.Q']
equiv_vars = ['source.P', 'source.Q']

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

identificator = greybox.DifferentialEvolution(model=equivmodel, measurements=yref,
                                     ref_vars=ref_vars, target_vars=equiv_vars,
                                     param_map=param_map)
if True:
    results = list(de(identificator.cost_fcn, [(uub, llb) for uub, llb in zip(lb,ub)]))
    print(results[-1])
    
# Setting the weights
weights = 0.01*np.ones(len(yref))
idxs = np.logical_and(yref.index > 3.28,yref.index < 8)
weights[idxs] = 1
identificator.weights = np.append(weights, weights)

bounds = [(lbv, ubv) for lbv, ubv in zip(lb,ub)]
ident_opts = {'x0':None, 'disp':True, 'bounds':bounds, 'updating':'deferred'}

###################### RUNNING THE IDENTIFICATION ######################
results = identificator.identify(**ident_opts)
theta_hat = np.array(results.x)

########################### PRINTING RESULTS ###########################

print(f'======= Estimated parameters  =======\n\n')
print(''.join([f'{name} = {theta_hat[idx]}\n' for idx, name in param_map.items()]))
print(f'=====================================')
print(f'Number of cost function calls: {identificator.i}\n\n')