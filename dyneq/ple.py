import numpy as np
import copy
from scipy.stats import chi2 as chi2
import dyneq

class ParPLEstimator(dyneq.identification.IdentificationProblem):
    def __init__(self, pplexperiment, model, experiments=[]):
        super().__init__(experiments, model)
        self.pplexperiment = pplexperiment

    def solve(self, indices, opt):
        if len(indices) > 1:
            if not opt.multicore:
                results = {}
                for idx in tqdm(indices):
                    results[idx] = self.estimate_pl(idx, opt)
                return results
            else:
                pass
                # NEEED TO FINISH THIS

    def estimate_pl(self, idx, opt):
        res = PLResults(opt.param_map[idx], idx)

        self._evaluate_branch('positive', idx, res, opt)
        self._evaluate_branch('negative', idx, res, opt)
        
        res.sort()
        return res

    def _evaluate_branch(self, side, theta_idx, res, opt):
        if side == 'positive':
            sign = 1
        elif side == 'negative':
            sign = -1
        else:
            raise ValueError('Incorrect PLE side option chosen.')
        pass
        # Initialization
        theta_prev = opt.theta_hat # Start from identified parameters
        dtheta_idx = sign*opt.max_step[theta_idx]/10 # Start with this change
        # Adjust the parameter map
        param_map = _adjust_map(theta_idx, opt.param_map)
        self.identificator.param_map = param_map
        
        for _ in range(opt.sample_size[theta_idx]):

            # Adjust the change
            self.identificator.param_map = opt.param_map
            dtheta, dtheta_idx = self.step_adjust(theta_idx, theta_prev, dtheta_idx)
            self.identificator.param_map = param_map
            print('Evaluating dtheta_idx: ', dtheta_idx)
            
            # IF not able to compute new theta end
            if np.any(np.isnan(dtheta)):
                break
            
            # Compute theta with the adjusted parameter step
            theta_trial = theta_prev + dtheta
            
            # Update the parameters in the model
            self.identificator.model.set_params(theta_trial, opt.param_map)
            
            identres = self.identificator.identify(x0=theta_prev[np.arange(len(theta_prev))!=theta_idx], jac='3-point', method='trf') # CHECK WHICH PARAMETER AND BOUNDS
            theta_prev = np.insert(identres.x, theta_idx, theta_trial[theta_idx])
            
            self.identificator.param_map = opt.param_map
            PL_theta = self.cost_fcn(theta_prev)
            res.PL = np.append(res.PL, PL_theta)
            res.theta = np.append(res.theta, theta_prev[theta_idx])
            
            if (self.cost_fcn(theta_prev) - opt.hat_cost_fcn) > opt.dchi2:
                break
            pass
        
        # Restore the parameter map
        self.identificator.param_map = opt.param_map
        
        return res      

class PLEstimator:
    
    def __init__(self, model, cost_fcn, identificator, opt):
        self.model = model
        self.cost_fcn = cost_fcn
        self.identificator = identificator
        self.opt = opt
        pass
         
    
    def estimate(self, theta_idx):
        '''
        Estimate a Profile Likelihood for a given parameter.
        '''
        # Creating a result object
        res = PLResults(self.opt.param_map[theta_idx], theta_idx)
        
        # Evaluating PL to the negative and positive of the optimized parameter
        self._evaluate_branch('positive', theta_idx, res)
        self._evaluate_branch('negative', theta_idx, res)
        
        # Add the optimum values
        res.theta = np.append(res.theta, self.opt.theta_hat[theta_idx])
        res.PL = np.append(res.PL, self.opt.hat_cost_fcn)
        
        # Sort the results
        res.sort()
        
        return res
    
    def _evaluate_branch(self, side, theta_idx, res):
        '''
        Performs one-sided evaluation from $\hat{\theta}$.
        
        Parameters
        ----------
        side : ['positive', 'negative'] 
            Defines 
        '''
        if side == 'positive':
            sign = 1
        elif side == 'negative':
            sign = -1
        else:
            raise ValueError('Incorrect PLE side option chosen.')
            
        # Initialization
        theta_prev = self.opt.theta_hat # Start from identified parameters
        dtheta_idx = sign*self.opt.max_step[theta_idx]/10 # Start with this change
        # Adjust the parameter map
        param_map = _adjust_map(theta_idx, self.opt.param_map)
        self.identificator.param_map = param_map
        
        for _ in range(self.opt.sample_size[theta_idx]):

            # Adjust the change
            self.identificator.param_map = self.opt.param_map
            dtheta, dtheta_idx = self.step_adjust(theta_idx, theta_prev, dtheta_idx)
            self.identificator.param_map = param_map
            print('Evaluating dtheta_idx: ', dtheta_idx)
            
            # IF not able to compute new theta end
            if np.any(np.isnan(dtheta)):
                break
            
            # Compute theta with the adjusted parameter step
            theta_trial = theta_prev + dtheta
            
            # Update the parameters in the model
            self.identificator.model.set_params(theta_trial, self.opt.param_map)
            
            identres = self.identificator.identify(x0=theta_prev[np.arange(len(theta_prev))!=theta_idx], jac='3-point', method='trf') # CHECK WHICH PARAMETER AND BOUNDS
            theta_prev = np.insert(identres.x, theta_idx, theta_trial[theta_idx])
            
            self.identificator.param_map = self.opt.param_map
            PL_theta = self.cost_fcn(theta_prev)
            res.PL = np.append(res.PL, PL_theta)
            res.theta = np.append(res.theta, theta_prev[theta_idx])
            
            if (self.cost_fcn(theta_prev) - self.opt.hat_cost_fcn) > self.opt.dchi2:
                break
            pass
        
        # Restore the parameter map
        self.identificator.param_map = self.opt.param_map
        
        return res      
    
    def step_adjust(self, idx, theta_prev, dtheta_idx_prev):
        '''
        Adjust the step in PL estimation. 
        
        Parameters
        ----------
        idx
            Index of the parameter whose PL is being estimated.
        theta_prev
            The previous set of parameters.
            
        
        '''
        # Some presetting
        step_factor = self.opt.step_factor[idx]
        dchi2 = self.opt.dchi2
        
        chi2last = self.cost_fcn(theta_prev)
        # Defining a vector of changes 
        dtheta = np.zeros(theta_prev.shape) 
        dtheta[idx]= dtheta_idx_prev
        
        # Making a copy for some reason ?
        dtheta_idx_new = copy.deepcopy(dtheta_idx_prev)
                
        while True:
            # Check if the parameter hit the limts
            hit_ub = theta_prev[idx] + dtheta[idx] >= self.opt.ub[idx]-self.opt.min_step[idx] # Parameter hit upper bound
            hit_lb = theta_prev[idx] + dtheta[idx] <= self.opt.lb[idx]-self.opt.min_step[idx] # Parameter hit lower bound

            if hit_ub or hit_lb:
                print(f'Param hit bound at {dtheta_idx_new}')
                # If it's hit then reduce it by a factor and
                dtheta_idx_new /= step_factor
                # If the reduced factor is smaller by the minimum step size then end
                if abs(dtheta_idx_new) < self.opt.min_step[idx]:
                    if dtheta_idx_prev > 0:
                        lbub = 'upper'
                    else:
                        lbub = 'lower'
                    dtheta = np.empty(self.opt.theta_hat.shape)
                    dtheta.fill(np.nan)
                    return dtheta, np.nan
            else:
                # Evaluate the cost function with new parameters
                chi2trial = self.cost_fcn(theta_prev+dtheta)
                if (chi2trial - chi2last > dchi2*self.opt.relchi2stepincrease[idx]):        
                    dtheta_idx_new /= step_factor
                    if np.abs(dtheta_idx_new) < self.opt.min_step[idx]:
                        return dtheta, dtheta_idx_new * step_factor
                    return dtheta, dtheta_idx_new
                else:
                    if ((chi2trial - chi2last) <= dchi2*self.opt.relchi2stepincrease[idx]) and (chi2trial-chi2last) > 0:
                        dtheta_idx_new *= step_factor
                        if np.abs(dtheta_idx_new) > self.opt.max_step[idx]:
                            dtheta_idx_new = self.opt.max_step[idx] * np.sign(dtheta_idx_new)
                    return dtheta, dtheta_idx_new
            dtheta[idx] = dtheta_idx_new

class PLResults:
    
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx
        self.PL = np.array([])
        self.theta = np.array([])
    
    def sort(self):
        sorted_idxs = self.theta.argsort()
        self.theta = self.theta[sorted_idxs]
        self.PL = self.PL[sorted_idxs]       
    
    def plot(self):
        pass
        
class PLOptions:
    '''
    Attributes
    ----------
    theta_hat : list
        Originally identified set of parameters.
    sample_size : list
        One-sided sample size of each profile likelihood.
    '''
    def __init__(self, theta_hat, hat_cost_fcn, param_map, ub, lb):
        self.theta_hat = theta_hat
        self.hat_cost_fcn = hat_cost_fcn
        self.param_map = param_map
        self.ub = ub
        self.lb = lb
        self.sample_size = 100 * np.ones(theta_hat.shape, dtype=int)
        self.relchi2stepincrease = 0.1 * np.ones(theta_hat.shape)
        self.max_step = 2*(ub-lb)/self.sample_size
        self.min_step = np.ones(theta_hat.shape)*5*1e-4
        self.step_factor = 1.5*np.ones(theta_hat.shape)
        self.alpha_level = 0.02
        self.dchi2 = chi2.ppf(1-self.alpha_level, 1)


def _adjust_map(idx, param_map):
    '''
    Returns an adjusted parameter map where an item idx is removed and
    all other items with index > idx are adjusted with -1.
    '''
    adjusted_map = {}
    for idx_old, param_name in param_map.items():
       if idx_old < idx:
           adjusted_map[idx_old] = param_name
       elif idx_old > idx:
           adjusted_map[idx_old-1] = param_name
    return adjusted_map

def cost_heatmap2d(cost_fcn, opt, theta, idx, idy):
    x, y = np.meshgrid(np.linspace(opt.lb[idx], opt.ub[idx], opt.sample_size[idx]),
                       np.linspace(opt.lb[idy], opt.ub[idy], opt.sample_size[idy]))
    z = np.empty(x.shape)
    for i in range(x.shape[0]):
        print(f'Doing row {i}')
        for j in range(x.shape[1]):
            theta_temp = copy.deepcopy(theta)
            theta_temp[idx] = x[i,j]
            theta_temp[idy] = y[i,j]
            z[i,j] = cost_fcn(theta_temp)
    return x, y, z