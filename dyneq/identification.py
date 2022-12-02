import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt


class IdentificationProblem:
    '''
    Consists of experiments, runs experiments and creates residuals or cost functions?
    '''
    def __init__(self, experiments, model):
        self.experiments = experiments
        self.model = model
    
    def residuals(self, theta):
        return np.concatenate([experiment.run(self.model, theta) for experiment in self.experiments])

    def normalized_residuals(self, theta):
        resid = self.residuals(theta)
        return resid/np.sqrt(len(resid))

    def mse(self, theta):
        nresid = self.normalized_residuals(theta)
        return np.dot(nresid, nresid)

    def se(self, theta):
        resid = self.residuals(theta)
        return np.dot(resid, resid)

    def solve(self, algo, ident_opt):
        if algo == 'scipy_least_squares':
            return op.least_squares(self.normalized_residuals, **ident_opt)
        elif algo == 'scipy_minimize':
            ident_opt['x0'] = self._fix_initial_guess(ident_opt['x0'], ident_opt['bounds'])
            return op.minimize(fun=self.mse, **ident_opt)        
        elif algo == 'scipy_de':
            return op.differential_evolution(self.mse, **ident_opt)
        elif algo == 'scipy_basin':
            return op.basinhopping(self.mse, **ident_opt)

    @staticmethod
    def _fix_initial_guess(x0, bounds):
        for idx, x in enumerate(x0):
            if x < bounds[idx][0]:
                x0[idx] = bounds[idx][0]
            elif x>bounds[idx][1]:
                x0[idx] = bounds[idx][1]
        return x0

class PPLIdentificationProblem(IdentificationProblem):
    
    def __init__(self, pplexperiment, experiments, model):
        super().__init__(experiments, model)
        self.pplexperiment = pplexperiment
        self.z = None

    def residuals(self, theta):
        if self.z is None:
            raise Exception('Cannot run the experiments because variable z is not set!')
        return np.append(super().residuals(theta), self.pplexperiment.run(self.model, theta, self.z))

class Experiment:
    '''
    Sets the environment to run a provided model. Given inputs, reference measurements, weights,
    and variable names to compare it computes weighted residuals.
    '''
    def __init__(self, yref, inputs, reference, target, simulation_settings=None, weights = None,
                starting_param = None, starting_param_map = None):
        self.yref = yref
        self.inputs = inputs
        self.reference = reference
        self.target = target
        self.weights = weights
        self.simulation_settings = simulation_settings
        self.starting_param = starting_param
        self.starting_param_map = starting_param_map
    
    def prepare_sim(self, model):
        model.set_params(self.starting_param, self.starting_param_map)

    def run(self, model, theta):
        y = self._execute_simulation(model, theta)
        return self._residuals(y)
    
    def _execute_simulation(self, model, theta):
        if self.starting_param is not None:
            self.prepare_sim(model)
        model.inputs = self.inputs
        model.set_params(theta)
        y = model.simulate(self.simulation_settings)
        if len(y) < 0 or y.isnull().any().sum()>0 or y.isna().any().sum()>0:
            return None
        return y
    
    def _residuals(self, y):
        # Allocate memory for the residuals
        ncols = len(self.reference)
        nrows = len(self.yref)
        # If the simulation failed return residuals with high value
        if y is None:
            return 1e20*np.ones(ncols*nrows)
        resid = np.empty(ncols*nrows)
        # Compute the regular residuals
        for col, variable_names in enumerate(zip(self.reference, self.target)):
            resid[col*nrows:(col+1)*nrows] = self.yref[variable_names[0]].values - y[variable_names[1]].values
        if self.weights is not None:
            resid *= self.weights
        return resid

class PPLExperiment(Experiment):
    
    def __init__(self, pl_idx, pl_name, sigma,  inputs, reference, target, 
                simulation_settings, yref = None, weights = None, include_all_residuals=False):
        super().__init__(yref, inputs, reference, target, simulation_settings, weights)
        self.pl_idx = pl_idx
        self.pl_name = pl_name
        self.sigma = sigma
        self.z = None
        self.include_all_residuals = include_all_residuals
    
    def run(self, model, theta):
        y = self._execute_simulation(model, theta)
        if y is not None:
            self.zpred_new = y[self.pl_name].values[self.pl_idx]
        else:
            self.zpred_new = 1e20
        if self.include_all_residuals:
            if self.z is None:
                return np.append(self._residuals(y), 0)
            else:
                return np.append(self._residuals(y), (self.z-self.zpred_new)/self.sigma)
        else:
            if self.z is None:
                return np.array([0])
            else:
                return np.array([(self.z-self.zpred_new)/self.sigma])

def compute_fit(y, yref):
    return np.sqrt(1/len(y)*np.linalg.norm(y-yref)**2)
    #return 1 - np.linalg.norm(yref-y)/np.linalg.norm(yref-np.mean(yref))

def compute_rmse(y, yref):
    return np.sqrt(1/len(y)*np.linalg.norm(y-yref)**2)