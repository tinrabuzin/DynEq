import numpy as np
import os
from scipy.stats import chi2 as chi2
from scipy import interpolate
import multiprocessing as mp
import pickle
from tqdm import tqdm
import pandas as pd
import dyneq.utils.draw as draw
import dyneq.identification
import scipy.optimize as op

class PPLEstimationProblem(dyneq.identification.IdentificationProblem):
    def __init__(self, pplexperiment, model, experiments=[]):
        super().__init__(experiments, model)
        self.pplexperiment = pplexperiment

    def solve(self, var_names, indices, opt):
        if len(var_names) > 1:
            if not opt.multicore:
                results = {}
                for var_name, idx in tqdm(zip(var_names, indices)):
                        results[idx] = self.estimate_ppl(var_name, idx, opt)
                return results    
            else:
                with mp.Pool(opt.nproc) as p:
                    arg_parall = [(var_name, idx, opt) for var_name, idx in zip(var_names, indices)]
                    return list(tqdm(p.imap(self._unpack, arg_parall), total=len(indices)))
        elif len(indices) == 1:
            return {indices[0]:self.estimate_ppl(var_names[0], indices[0], opt)}
    
    def _unpack(self, args):
        return self.estimate_ppl(*args)

    def estimate_ppl(self, var_name, idx, opt):
        # Create the objects where results will be stored for the current point
        res = PLResults(var_name, idx)
        # Define the point and sigma for which PPL will be estimated
        self.pplexperiment.pl_idx = idx
        self.pplexperiment.pl_name = var_name
        self.pplexperiment.sigma = opt.sigma
        # Run the first identification with the validation point
        # and prepare options for PPL estimation
        hat_res = self.identify(opt.solver, opt.solver_opts)
        opt.theta_hat = hat_res.x
        opt.solver_opts['x0'] = hat_res.x
        opt.hat_cost_fcn = self.se(opt.theta_hat)
        self.zpred_orig = self.pplexperiment.zpred_new
        
        self._evaluate_branch('negative', res, opt)
        self._evaluate_branch('positive', res, opt)
        # Sort the results
        res.sort()
        # Save the results
        if opt.save_intermediate:
            res_path = os.path.join(opt.results_path, f'pl_{self.model.model_name}_{var_name}_{idx}.p')
            pickle.dump(res, open(res_path ,'wb'))
        return res

    def identify(self, solver, solver_opts):
        return super().solve(solver, solver_opts)

    def residuals(self, theta):
        if len(self.experiments) > 0:
            residuals = [experiment.run(self.model, theta) for experiment in self.experiments]
        else:
            residuals = []
        residuals.append(self.pplexperiment.run(self.model, theta))
        return np.concatenate(residuals)

    def _evaluate_branch(self, side, res, opt):
        if side == 'positive':
            sign = 1
        elif side == 'negative':
            sign = -1
        else:
            raise ValueError('Incorrect PLE side option chosen.')

        chi2_new = 0        
        first_step = True
        for i in range(opt.maxsteps):
            # Saving the values from the previous steps or defining the
            # values if the step is first.
            if first_step:
                theta = opt.theta_hat
                dz = opt.firststep 
                z_old = self.zpred_orig
                chi2_old = 0
                first_step = False # Next step is not first anymore
            else:
                z_old = self.pplexperiment.z
                # Define the new step to take
                if opt.stepping_method == 'past_based':
                    dz = self.pastbased_step_adjustment(merit_diff, dz, opt)
                elif opt.stepping_method == 'dynamic':
                    dz = self.dynamic_step_adjustment(i, sign, z_old, merit_diff, dz, self.pplexperiment.zpred_new, opt)
            self.pplexperiment.z = z_old + sign*dz
            # Identify the new set of parameters with z_new
            opt.solver_opts['x0'] = theta
            ident_res = self.identify(opt.solver, opt.solver_opts)
            theta_old = theta
            theta = ident_res.x
            cost_new = self.se(theta)
            chi2_new = cost_new-opt.hat_cost_fcn
            ppl = chi2_new - ((self.pplexperiment.zpred_new-self.pplexperiment.z)/opt.sigma)**2
            if opt.stop_criteria == 'ppl':
                merit_diff = ppl - chi2_old
                chi2_old = ppl
            else:
                merit_diff = chi2_new - chi2_old
                chi2_old = chi2_new

            if chi2_new < 0:
                print(f'Not in global optimum -  {chi2_new}')
                print(ident_res.x)
                with open(f'theta_{res.idx}.p', 'wb') as file:
                    pickle.dump(ident_res.x, file)
            # Compute the point in prediction profile likelihood
            try:
                if np.abs(res.zhat[-1] - self.pplexperiment.zpred_new) < 1e-16:
                    continue
            except IndexError:
                pass

            # Store the results
            res.z = np.append(res.z, self.pplexperiment.z)
            res.vpl = np.append(res.vpl, chi2_new)
            res.ppl = np.append(res.ppl, ppl)
            res.zhat = np.append(res.zhat, self.pplexperiment.zpred_new)
            #res.th = np.append(res.th, ident_res)
            #pickle.dump(res, open(os.path.join(opt.results_path, f'pl_{self.model.model_name}_{self.pplexperiment.pl_name}_{self.pplexperiment.pl_idx}.p'),'wb'))
            print(f'Now at {res.idx} with PPL:\n {res.ppl}')
            if opt.stop_criteria == 'vpl':
                if chi2_new >= opt.chi2max:
                    break
            elif opt.stop_criteria == 'ppl':
                if res.ppl[-1] >= opt.chi2max:
                    break
            if i > opt.maxsteps:
                print(f'Exiting early for idx: {res.idx}')
                break
    
    @staticmethod
    def pastbased_step_adjustment(chi2_dif, dz, opt):
        '''Adjust the step in PL estimation.  
        '''
        maxstepsize       = opt.maxstepsize
        minstepsize       = opt.minstepsize
        stepfactor        = opt.stepfactor
        chi2dif_max       = opt.chi2dif_max
        chi2dif_min       = opt.chi2dif_min
        
        if abs(chi2_dif) > chi2dif_max:
            if dz/stepfactor > minstepsize:
                dz /= stepfactor
            return dz
        elif abs(chi2_dif) < chi2dif_min:
            if dz*stepfactor < maxstepsize:
                dz *= stepfactor
            return dz
        else:
            return dz

    def dynamic_step_adjustment(self, i, sign, z_old, chi2_dif, dz, pred_old, opt):
        chi2dif_max = opt.chi2dif_max
        maxstepsize = opt.maxstepsize
        minstepsize = opt.minstepsize
        stepfactor  = opt.stepfactor
        sigma       = opt.sigma
        
        # Theoretical upper limit for chi2 change
        chi2limit = self._chi2limit(sign, dz, sigma, z_old, pred_old)
        
        if chi2limit > chi2dif_max:
            while chi2limit > chi2dif_max:
                if dz/stepfactor > minstepsize:
                    dz /= stepfactor
                else:
                    dz = minstepsize
                    break
                chi2limit = self._chi2limit(sign, dz, sigma, z_old, pred_old)
        else:
            while chi2limit <= chi2dif_max:
                if dz*stepfactor < maxstepsize:
                    dz *= stepfactor
                else:
                    dz = maxstepsize
                    dz *= stepfactor
                    break
                chi2limit = self._chi2limit(sign, dz, sigma, z_old, pred_old)
            dz /= stepfactor

        if not ((dz == maxstepsize) or (dz == minstepsize)) and (abs(chi2_dif) < opt.chi2dif_min) and i != 1:
               dz = stepfactor*dz
               
        return dz

    @staticmethod
    def _chi2limit(sign, dz, sigma, z_old, pred_old):
        return (sign*dz/sigma**2)*(2*(z_old-pred_old)+sign*dz)

class PLResults:
    """A class that stores the VPL and PPL estimation results.

    Attributes
    ----------
    name : str
        Name of the variable for which the estimation is to be done.
    idx : int
        Index at which the estimation is to be done. 
    z
        Points at which cost function was estimated.
    vpl
        The value of the cost function at points z, i.e. VPL.
    ppl
        Prediction profile likelihood at z.
    """
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx
        self.vpl = np.array([])
        self.z = np.array([])
        self.ppl = np.array([])
        self.zhat = np.array([])
        self.th = np.array([])
        self.y = pd.DataFrame()
        self.ti = -1

    def sort(self):
        """Sorts the estimation results with z as x-axis.
        """        
        sorted_idxs = self.z.argsort()
        self.z = self.z[sorted_idxs]
        self.vpl = self.vpl[sorted_idxs]
        #self.th = self.th[sorted_idxs]
        sorted_idxs = self.zhat.argsort()
        self.zhat = self.zhat[sorted_idxs]
        self.ppl = self.ppl[sorted_idxs]

    def plot(self, option='ppl', axs = None):
        return draw.plot_plres(self, self.ti, option='ppl', axs=axs)

class PLOptions:
    '''
    Attributes
    ----------
    theta_hat : list
        Originally identified set of parameters.
    alpha_level : float
        Significance level 
    chi2max : float
        Maximum Chi2 profile value
    maxsteps : int
    Maximum number of steps in the evaluation
    
    '''
    def __init__(self, theta_init, sigma, alpha_level, multicore=False, 
                nproc = 0, stop_criteria='ppl', save_intermediate=False, 
                solver='scipy_minimize', solver_opts={}, results_path='', stepping_method='past_based'):
        self.theta_init = theta_init
        self.sigma = sigma
        self.chi2dif_max = 0.3
        self.chi2dif_min = 0.05
        self.alpha_level = alpha_level
        self.chi2max = chi2.ppf(1-self.alpha_level, 1)
        self.maxsteps = 100
        self.stepfactor = 1.5 # Step size adaption factor
        self.maxrange = 50*sigma # Maximal absolute range of profile in each direction
        self.maxstepsize = 5*sigma # Maximal absolute value of step
        self.minstepsize = 5*sigma/1000 # Minimal absolute value of step
        self.firststep = sigma/5 # First step from the optimum
        self.multicore = multicore
        self.save_intermediate = save_intermediate
        self.results_path = results_path
        self.solver = solver
        self.solver_opts = solver_opts
        if nproc == 0 and multicore:
            self.nproc = mp.cpu_count()
        else:
            self.nproc = nproc
        if (stop_criteria != 'ppl') and (stop_criteria != 'vpl'):
            raise ValueError('Select one of the valid stop criteria [ppl/vpl]!')
        else:
            self.stop_criteria = stop_criteria
        if (stepping_method != 'past_based') and (stepping_method != 'dynamic'):
            raise ValueError('Select one of the valid stop criteria [past_based/dynamic]!')
        else:
            self.stepping_method = stepping_method

def extract_ci(pl_results, pl_option='vpl'):
    alpha_level = 0.05
    chi2max = chi2.ppf(1-alpha_level, 1)

    z_max = np.zeros(len(pl_results))
    z_min = np.zeros(len(pl_results))
    t_z = np.empty(len(pl_results))

    for idx, pl in enumerate(pl_results):
        ppl_interp = interpolate.interp1d(pl.zhat, pl.ppl-chi2max,'linear', fill_value='extrapolate')
        z_min[idx], z_max[idx] = op.newton(ppl_interp, pl.zhat[1]), op.newton(ppl_interp, pl.zhat[-2])
        t_z[idx] = pl.ti
    
    return pd.DataFrame(data={'zmax':z_max, 'zmin':z_min}, index=t_z)
