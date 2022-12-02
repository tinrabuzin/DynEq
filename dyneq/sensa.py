import numpy as np
import multiprocessing as mp

    
def execute(model, param_map, param_range, multicore=True):
    '''
    Executes sensitivity analysis.
    '''
    if multicore:
        return _simulate_multicore(model, param_map, param_range)
    else:
        return _simulate_sequential(model, param_map, param_range)

def _simulate_multicore(model, param_map, param_range):
    arg_parall = [(model, param_map, param_values) for param_values in param_range]
    nprocs = mp.cpu_count()
    pool = mp.Pool(processes=nprocs)
    return pool.starmap(_simulate, arg_parall)
    
def _simulate_sequential(model, param_map, param_range):
    for param_values in param_range:
        _simulate(model, param_map, param_values)

def _simulate(model, param_map, param_values):
    model.set_params(param_values)
    return model.simulate()