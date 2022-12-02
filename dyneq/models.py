import os
import csv
import re
import copy
import warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import subprocess
from scipy.linalg import inv, expm, eig, eigvals, logm
import tempfile
import uuid
import atexit

try:
    import dyneq.pf as pf
    from dyneq.pf import pfhandle
except ImportError:
    warnings.warn('Failed importing PowerFactory dependencies!')
import dyneq.wavegen as wavegen
try:
   from OMPython import OMCSessionZMQ, ModelicaSystem
   import DyMat
   omchandle = None
except:
   warnings.warn('Failed importing OpenModelica dependencies!')
try:
   import fmpy
except ImportError:
   warnings.warn('Failed importing FMU dependencies!')

class AbstractModel:

    def __init__(self):
        raise NotImplementedError

    def simulate(self, t_start, t_stop, dt=None):
        raise NotImplementedError

    def linearize(self):
        raise NotImplementedError

    def set_params(self, *args):
        raise NotImplementedError

    def get_params(self, *args):
        raise NotImplementedError

class FMU(AbstractModel):

    def __init__(self, build_fmu, **kwargs):
        if build_fmu:
            self.fmu_path = self.build_fmu(**kwargs)
        else:
            self.fmu_path = kwargs['fmu_path']
        self.output = kwargs.pop('output', None)
        print('Creating a model with the following FMU:')
        fmpy.dump(self.fmu_path)
        self.model_description = fmpy.read_model_description(self.fmu_path)#, validate=False)
        self.unzipdir = fmpy.extract(self.fmu_path)
        self.fmu=fmpy.instantiate_fmu(self.unzipdir, self.model_description, 'ModelExchange')
        self.input_specs = kwargs.pop('input_specs', None)
        self._create_input(self.input_specs)
        self.postprocess = kwargs.pop('postprocess', None)
        self.start_values = {}


    def _create_input(self, input_specs=None):
        if input_specs is None:
            self.input = None
        else:
            with open(input_specs, newline='\n') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                names = next(reader)

            signal = np.genfromtxt(input_specs, delimiter=',', skip_header=1, dtype=[(name, np.double) for name in names])
            signal = signal[list(signal.dtype.names[:-1])]
            self.input = signal
            
        
    def simulate(self, experiment_settings={}):
        self.fmu.reset()
        results = fmpy.simulation.simulate_fmu(self.unzipdir, fmu_instance=self.fmu, model_description=self.model_description,
                                           validate=False, input=self.input, start_values=self.start_values, output=self.output, **experiment_settings)
        
        if self.postprocess is not None:
            return self.postprocess(results)
        else:
            return results
    
    def set_params(self, theta):
        if len(args) == 1:
            self.start_values = args[0]    
        elif len(args) == 2:
            self.start_values = {param_name:args[0][param_idx] for param_idx, param_name in args[1].items()}
        else:
            raise ValueError('Method parameters are not correct.')

        
    def build_fmu(self, **kwargs):
        if not omchandle:
            start_openmodelica()
        model_name = kwargs['model_name']
        model_file_path = kwargs['model_file_path']
        try:
            dependencies = kwargs['dependencies']
        except KeyError:
            dependencies = []
        # Change to target C++
        omchandle.sendExpression(f'setCommandLineOptions("--simCodeTarget=Cpp")')
        omchandle.sendExpression(f'setCommandLineOptions("--target=gcc")')
        omchandle.sendExpression(f'setCompiler("gcc")')
        omchandle.sendExpression(f'setCXXCompiler("g++")')
        # Change to where the FMU will be built
        omchandle.sendExpression(f'cd("{kwargs["build_path"]}")')
        # Use new instatiation method which helps with complex parameters
        omchandle.sendExpression('setCommandLineOptions("--replaceHomotopy=actual)')
        # Load Modelica library
        omchandle.sendExpression('loadModel(Modelica)')
        # Load all dependencies
        for dependency in dependencies:
            if not omchandle.sendExpression(f'loadFile("{dependency}")'):
                raise Exception('Could not load dependency!')
        # Load the model
        omchandle.sendExpression(f'loadFile("{model_file_path}")')
        build_cmd =  f'buildModelFMU({model_name},version="2.0",fmuType="me")'
        fmu_path = omchandle.sendExpression(build_cmd)
        if not fmu_path:
            raise Exception('Could not build FMU!')
        return fmu_path
    
    def set_params(self, *args):
        if len(args) == 1:
            params = args[0]    
        elif len(args) == 2:
            params = {param_name:args[0][param_idx] for param_idx, param_name in args[1].items()}
        else:
            raise ValueError('Method parameters are not correct.')
        self.model.parameters = params
    
class OMModel(AbstractModel):
    
    def __init__(self, existing, **kwargs):
        self.model_name = kwargs['model_name']
        if not omchandle and not existing:
            start_openmodelica()
        self._build_new(**kwargs)
        self.postprocess = kwargs.pop('postprocess', None)
        self._params = {}        
    
    def _build_new(self, **kwargs):
        '''Compiles the binary used for simulations.
        '''
        self.model_name = kwargs['model_name']
        self.model_file_path = kwargs['model_file_path']
        try:
            self.dependencies = kwargs['dependencies']
        except KeyError:
            self.dependencies = []
        try:
            self.input_specs = kwargs['input_specs']
        except KeyError:
            self.input_specs = ''
        self.output_specs = kwargs['output_specs']
        try:
            self.param_map = kwargs['param_map']
        except KeyError:
            self.param_map = None
        self._experiment_settings = kwargs['experiment_settings']
        # Change to where the binary will be built
        self.tempfolder = tempfile.TemporaryDirectory(prefix='omc_')
        self.binary_folder = self.tempfolder.name
        self.binary_path = os.path.join(self.binary_folder, self.model_name)
        omchandle.sendExpression(f'cd("{self.binary_folder}")')
        # Use new instatiation method which helps with complex parameters
        omchandle.sendExpression('setCommandLineOptions("--replaceHomotopy=actual")')
        # Load Modelica library
        omchandle.sendExpression('loadModel(Modelica)')
        # Load all dependencies
        for dependency in self.dependencies:
            if not omchandle.sendExpression(f'loadFile("{dependency}")'):
                raise Exception('Could not load dependency!')
        # Load the model
        omchandle.sendExpression(f'loadFile("{self.model_file_path}")')
        build_cmd = f'buildModel({self.model_name}, outputFormat="mat",'
        build_cmd += ','.join([f'{setting}={value}' for setting, value in kwargs['experiment_settings'].items() if setting != 'solver'])
        if 'solver' in kwargs['experiment_settings']:
            build_cmd += ','+f'method="{kwargs["experiment_settings"]["solver"]}"'+')'
        else:
            build_cmd += ')'
        # Build the model                 
        omchandle.sendExpression(build_cmd)

    def simulate(self, experiment_settings=None):
        '''
        Simulates the model by running the compiled binary.

        Notes
        -----
        The outputs are processed to obtained equidistant simulation points. This should be removed
        when https://trac.openmodelica.org/OpenModelica/ticket/5709 is fixed.
        '''
        command = self._simulation_command(experiment_settings)
         #
        try:
            subprocess.run(command,stdout=subprocess.DEVNULL,cwd=self.binary_folder, shell=True, timeout=2) #  
            mat = DyMat.DyMatFile(self.results_path)
        except:
            return pd.DataFrame(
            data=[[100,100]],
            columns=self.output_specs,
            index = [0]
        )
        os.remove(self.results_path)
        # Equivdistancing
        if experiment_settings is None:
            experiment_settings = self._experiment_settings
        if 'numberOfIntervals' in experiment_settings:
            newidx = np.round(np.linspace(experiment_settings['startTime'],experiment_settings['stopTime'],int(experiment_settings['numberOfIntervals'])),2)
        else:
            step = experiment_settings['stepSize']
            newidx = np.round(np.arange(experiment_settings['startTime'],experiment_settings['stopTime']+step, step), 2)
        oldidx = mat.abscissa(self.output_specs[0])[0]
        interp_data = np.empty((len(newidx),len(self.output_specs)))
        if len(oldidx) == 0:
            interp_data = 1e10*np.ones(interp_data.shape)
        else:    
            for colidx, var_name in enumerate(self.output_specs):
                interp_data[:, colidx] = np.interp(newidx, oldidx, mat.data(var_name))
        # Storing results in Pandas frame
        results = pd.DataFrame(
            data = interp_data,
            columns=self.output_specs,
            index=newidx
        )
        if self.postprocess is not None:
            return self.postprocess(results)
        else:
            return results
        
    def _simulation_command(self, experiment_settings):
        '''Creates the full simulation command used in the subprocess.
        '''
        self.results_path = os.path.join(self.binary_folder, f'results_{uuid.uuid4()}.mat').replace("\\","/")
        #self.results_path = os.path.join(self.binary_folder, f'results_{os.getpid()}.mat').replace("\\","/")
        flags_str = " ".join([f'-r={self.results_path}', 
                              self._override_flag(self.results_path, experiment_settings),
                              #'-noEventEmit',
                              '-lv=-stdout,-assert',
                              #'-newtonFTol=1e-14',
                              #'-newtonXTol=1e-14',
                              #'-nls=newton',
                              self._input_flag()])
        return f'{self.binary_path} {flags_str}'

    def _override_flag(self, results_path, experiment_settings):
        '''
        Creates an override flag to set output variable filter,
        experiment settings an the model's parameters.
        '''
        if self.output_specs:
            variable_filter = f'variableFilter={self._output_var_flag(self.output_specs)}'
        else:
            variable_filter = ''
            
        if experiment_settings:
            experiment_strings = [
            f'{str(variable)}={str(value)}' for variable, value in experiment_settings.items()]
        else:
            experiment_strings = [
            f'{str(variable)}={str(value)}' for variable, value in self._experiment_settings.items()    
            ]        ## FIX THIS CHECK
        if experiment_strings or variable_filter:
            override_flag = ",".join([variable_filter]+experiment_strings)
            for name,value in self._params.items():
                override_flag += f',{name}={value:.16f}'
            return f'-override={override_flag} '
        else:
            return ''    

    @property
    def inputs(self):
        return self.input_specs
    @inputs.setter
    def inputs(self, inputs):
        self.input_specs = inputs

    @staticmethod
    def _output_var_flag(output_variables):
        '''Creates REGEX for variable filtering in outputs.
        '''
        return "'|'".join(['^' + outvar.replace('.','\\.') + '$' for outvar in output_variables])

    def _input_flag(self):
        '''Creates a flag that is used to supply an input CSV file for the simulation command.
        '''
        if self.input_specs:
            return f'-csvInput={self.input_specs}'
        else:
            return ''

    def set_params(self, *args):
        '''Creates and stores a parameter dictionary which is set using a flag in the simulation
        '''
        if len(args) == 1:
            param_map = self.param_map
        else:
            param_map = args[1]
        new_params = {param_name:args[0][param_idx] for param_idx, param_name in param_map.items()}
        self._params.update(new_params)      
        
    def get_params(self, params):
        '''Returns the previously stored parameter dictionary
        '''
        return self._params
    
class PowerFactoryModel(AbstractModel):
    
    def __init__(self, project_name, project_folder=None, **kwargs):
        # Starting PowerFactory
        if pf.pfhandle is None:
            pf.start_powerfactory()
        pf.pfhandle.ResetCalculation() # Reseting any calculations from previous models
        # Activating the project
        self._project_name = project_name
        self.prj = pf.Project(project_name, project_folder)
        # Used to store the byproducts such as linearization results
        try:
            self._tempdir = kwargs['tempdir']
        except KeyError:
            self._tempdir = pf.pfhandle.GetSettings('tempdir')
        # Selecting which study case to use and storing it at self.stc
        try:
            self._stc_name = kwargs['study_case']
        except KeyError:
            pass
        else:
            self.prj.study_cases[self._stc_name].activate()
        self.stc = self.prj.study_cases.active_study_case
        if self.stc is None: ### CHECK THIS
            raise RuntimeError('No study case has been activated!')
        # Selecting which networks to be activated in the study case
        try:
            if type(kwargs['networks']) is list:
                self._act_nets = kwargs['networks'] # Store the selection
            else:
                self._act_nets = [kwargs['networks']]
            for name, net in self.prj.networks.items():
                if any(name == act_name for act_name in self._act_nets):
                    net.activate()
                else:
                    net.deactivate()
        except KeyError:
            pass
        # Defining the model outputs
        try:
            self._output_specs = kwargs['outputs']
        except KeyError:
            pass
        else:
            self.outputs = self._output_specs
        self.out_type = 'db' # Define if reading from ElmRes or ComRes
        # Defining the model inputs
        try:
            self._input_specs = kwargs['inputs']
        except KeyError:
            pass
        else:
            self.inputs = self._input_specs
    
    @staticmethod
    def set_params(*args):
        """ #### NOT FINISHED Maps a list of parameters into a dictionary.

        This is done so that it is possible to set the parameters of the
        equivalent method using :code:`pf.set_params()`.

        Parameters
        ----------
        param : list
            A list of equivalent model's parameters.
        param_map : dict
            A dictionary to map the parameter list. It is of the 
            following form
            :code:`{0:('element name 0', 'parameter name 0'),1:...}`

        Returns
        -------
        dict
            The dictionary of parameters of the following form:
            :code:`{'element name 0':{'parameter name 0':value, ... }',
            ...}`
        """
        if len(args) == 1:
            pf.set_params(args[0])
        elif len(args) == 2:
            theta = {param_specs[0]:{param_specs[1]:args[1][param_pos]} for param_pos, param_specs in args[0].items()}
            pf.set_params(theta)
        else:
            raise ValueError('Method parameters are not correct.')

    def setup(self):
        ### CHECK FOR PROJECT ACTIVATION
        if not self.prj.study_cases[self._stc_name] is self.prj.study_cases.active_study_case:
            self.prj.study_cases[self._stc_name].activate()
            self.stc = self.study_cases.active_study_case
        try:
            act_nets = self._act_nets
        except AttributeError:
            pass
        else:
            for name, net in self.prj.networks.items():
                if name in act_nets:
                    net.activate()
                else:
                    net.deactivate()
        try:
            self.outputs = self._output_specs
        except AttributeError:
            pass
        try:
            self.inputs = self._input_specs
        except AttributeError:
            pass

    def simulate(self, t_start, t_stop, dt, output=True):
        self.elmres.state = pf.ElmRes.NOTLOADED
        self.stc.simulate(t_start, dt, t_stop, 'rms')
        if output == True:
            return self.outputs

    @property
    def outputs(self):
        if self.out_type == 'db':
            return self.elmres.outputs
        elif self.out_type == 'file':
            return self.comres.read_results()
        else:
            raise Exception('Output type not specified!')
        
    @outputs.setter
    def outputs(self, specs):
        try:
            elmres_name = specs['ElmRes']
        except KeyError:
            elmres = self.elmres
        else:
            try:
                elmres = self.stc.elmres[elmres_name] ### FINISH UP
            except KeyError:
                elmres = self.stc.create_elmres(elmres_name)
            self.elmres = elmres
        try:
            self.elmres.naming = specs['naming']
        except KeyError:
            self.elmres.naming = 'reduced'
        self.elmres.variables = specs['variables']
        self.stc.inc.p_resvar = self.elmres.pfobject
        try:
            comres_name = specs['ComRes']
        except KeyError:
            pass
        else:
            try:
                comres = self.stc.comres[comres_name] ### FINISH UP
            except KeyError:
                comres = self.stc.create_comres(comres_name)
            self.comres = comres
            self.comres.define_outputs(specs['variables'], self.elmres, specs['filepath'])
        
    @property
    def inputs(self):
        return self._input_specs
    
    @inputs.setter
    def inputs(self, inputs):  #### FORGOT TO SET self._inputs_specs
        for input_specs in inputs:
            net = input_specs['network']
            try:
                # Defining ElmFile to be input
                elmfile_name = input_specs['ElmFile']
                elmfile = self.prj.networks[net].elmfiles[elmfile_name]
                if input_specs['source'] == 'generate':
                    wave_specs = input_specs['wave_specs']
                    signal = wavegen.wavetype[wave_specs['type']](**wave_specs)
                    try:
                        filepath = input_specs['filepath']
                    except KeyError:
                        filepath =  os.path.join(self._tempdir,
                                            f'signal_{elmfile_name}.csv')
                    elmfile.create_file(signal, filepath)
                elif input_specs['source'] == 'ext_file':
                    elmfile.f_name = input_specs['filepath']   
                else:
                    raise ValueError('Input specifications not valid.')
            except KeyError:
                # Defining IntMat to be input
                intmat_name = input_specs['IntMat'] 
                intmat = self.prj.networks[net].intmats[intmat_name]
                wave_specs = input_specs['wave_specs']
                intmat.signal = wavegen.wavetype[wave_specs['type']](**wave_specs)
                
    def linearize(self, inputs_c, inputs_bd, outputs, ref_state=''):
        """Returns the `A`, `B`, `C` and `D` matrices of a numerically 
        linearized system, i.e. the matrices of the following system.
        
        .. math:: \dot{x} = Ax + Bu
        .. math:: y = Cx + Du
        
        Compute a column of B assuming steady state after a disturbance:  
        .. math: 0 = \dot{x} = Ax + B_iu_i \Rightarrow B_i = -Axu^{-1}
        .. math: D_i = (y-Cx)u^{-1}

        Parameters
        ----------
        inputs_c : list
            A list containing dictionaries that specify all the input 
            waveforms for computation of the matrix `C`.
        inputs_bd : list
            A list containing dictionaries that specify all the input 
            waveforms for computation of matrices `B` and `D`.
        outputs: dict[str, list]
            Dictionary of outputs where key is the name of an object in 
            PowerFactory and the value is a list of the correspodning 
            signals taken as an output of the system.
        outputs_idx: dict[str, int]
            Dictionary where a key is of the following form - 
            :code:`'ObjectName.ObjectType\\SignalName`
        modalpath: str, optional
            A path to a folder where PowerFactory will store the modal
            analysis results.
        ref_state: str
            This can be used to eliminate a reference state which causes
            the A matrix to be nonsingular. The state has to be
            specified using the following form 
            :code:`'ObjectName.ObjectType\\SignalName`
 
            
        Notes
        -----
        First, the system is perturbed from the steady state, and the
        unforced system response is recorded (both the system states and
        the outputs). The matrix C is found by solving the following least-squares
        problem 
        
        .. math:: C = yX^T[XX^T]^{-1}
            
        Once this is solved, a step change is introduced to each of the 
        inputs separately, and the values at the end of the simulation 
        time are assumed to be at the steady state. The colums of the
        matrices B and D are then computed as
        
        .. math:: 0 = \dot{x} = Ax + B_iu_i \Rightarrow B_i &= -Axu_i^{-1}
        
                  D_i &= (y-Cx)u_i^{-1}
        
        """
        # Taking out the current study case
        stc = self.stc
        # Setting the inputs before modal analysis
        self.inputs = inputs_c
        # Running the PowerFactory modal analysis
        self.stc.modal_analysis(self._tempdir)

        # Getting the states dictionaries
        states = copy.deepcopy(self.states)
        
        # Getting the dimensions 
        n_x = len(states['states_idx']) # number of states 
        n_y = len(outputs['outputs_idx']) # number of outputs 
        n_u = len(inputs_c) # number of inputs 
        
        # Remove a state in the case it's, e.g. a reference angle
        if ref_state:
            ref_idx = states['states_idx'].pop(ref_state) # save index for later
            model_name, state_variable = re.match(r'(.+)\\(s:.+)',ref_state).groups()
            states['variables'][model_name].remove(state_variable)
            if not states['variables'][model_name]:
                states['variables'].pop(model_name)
        
        # Prepare the outputs used for linearization
        states_outputs = {'ElmRes' : 'ModalComp', 'variables':{}}
        states_outputs['variables'].update(states['variables'])
        states_outputs['variables'].update(outputs['variables'])
        temp_elmres = self.elmres
        self.outputs = states_outputs
        
        # Get a list of ordered states
        ordered_states = [state for state, _ in sorted(states['states_idx'].items(), key = lambda item: item[1])]
        # Get a list of  ordered outputs    
        ordered_outputs = [output for output, _ in sorted(outputs['outputs_idx'].items(), key = lambda item: item[1])]
        
        
        # Getting the A matrix from PowerFactory results
        A = self.a_mat
        if ref_state is not None: # Removing a referece state if needed
            A = np.delete(np.delete(A,ref_idx,0),ref_idx,1)
            n_x -= 1

        # Simulate a step change on all inputs to get the matrix C
        wave_specs = inputs_c[0]['wave_specs']
        self.stc.inc.iopt_fastout = 0
        self.stc.inc.dtout = wave_specs['step']
        res = self.simulate(wave_specs['tstart'], wave_specs['tstop'], wave_specs['step'])
        res -= res.iloc[0,:] # Centering with respec to steady-state
        X = res[ordered_states].loc[(res.index >0.005) & (res.index < 0.5) | (res.index >0.503)].values.T 
        y = res[ordered_outputs].loc[(res.index >0.005) & (res.index < 0.5) | (res.index >0.503)].values.T
        C = np.linalg.solve(X@X.T, X@y.T).T
        
        # Computing the matrices B and D
        # Initializing of the matrices B and D
        B = np.empty((n_x,n_u)) 
        D = np.empty((n_y,n_u))
        
        for enab_idx, enab_inpt in enumerate(inputs_bd):
            inpts_specs = [enab_inpt] ### CHECK
            for disabl_idx in (idx for idx in range(n_u) if idx != enab_idx):
                disabl_inpt = copy.deepcopy(inputs_bd[disabl_idx])
                disabl_inpt['wave_specs']['type'] = 'const'
                inpts_specs.append(disabl_inpt)
            self.inputs = inpts_specs
            wave_specs = enab_inpt['wave_specs']
            res = self.simulate(wave_specs['tstart'],wave_specs['tstop'],wave_specs['step'])
            res -= res.iloc[0,:] # Center with respect to steady-state
            end_state = res[ordered_states].values[-1,:]
            end_output = res[ordered_outputs].values[-1,:]
            B[:, enab_idx] = -A@end_state/(-wave_specs['deltay'])
            D[:,enab_idx] = (end_output-C@end_state)/(-wave_specs['deltay'])
        
        self.elmres = temp_elmres
        self.stc.inc.p_resvar = self.elmres.pfobject
        return A, B, C, D
    
    @property
    def states(self):
        if hasattr(self, '_states'):
            return self._states
        else:
            self._states = self.get_states()
            return self._states
        
    def get_states(self):
        """
        Returns two dictionaries containing states and their indices
        used by PowerFactory. They are generated from modal computation
        results so they need to be generated first.
                
        Parameters 
        ----------
        modalpath : str
            A path to the folder where PowerFactory modal results are 
            stored.
        
        Returns
        -------
        states : dict [str, list]
            States is a dictionary where states are grouped by a network 
            element.
        states_idx : dict [str, int]
            A dictionary where every state as a key has a corresponding 
            index in PowerFactory as a value.
        """
        states = {}
        states['variables'] = defaultdict(list)
        states['states_idx'] = {}
        # These two dictionaries are read from VariableToIdx_Amat.txt 
        # generated as an output of modal analysis by PowerFactory.
        modalpath = self.stc.mod.dirMatl
        with open(os.path.join(modalpath,'VariableToIdx_Amat.txt'), 'r') as fp:
            for line in fp:
                match = re.match(r'^\s*([0-9]+)\s+([\S\s]+.\S+)\s+"(\S+)"', line)
                try:
                    state_idx, model_name, state_variable = match.groups()
                except (AttributeError, ValueError):
                    pass
                else:
                    states['states_idx'][f'{model_name}\\s:{state_variable}'] = int(state_idx)-1
                    states['variables'][model_name].append(f's:{state_variable}')
        return states

    @property
    def a_mat(self):
        if hasattr(self, '_a_mat'):
            return self._a_mat
        else:
            self._a_mat = self.read_amat()
            return self._a_mat

    def read_amat(self):
        """Returns the matrix A read from the file Amat.mtl which is the 
        output form PowerFactory modal analysis.
    
        Returns
        -------
        ndarray
            Returns the A matrix of the linearized system.
        """
        n = len(self.states['states_idx'])
        Amat = np.zeros((n,n))
        with open(os.path.join(self.stc.mod.dirMatl,'Amat.mtl')) as fp:
            for line in fp:
                row, col, value = line.split()
                Amat[int(row)-1, int(col)-1] = float(value)
        return Amat

    def linearize2(self, inputs, outputs, ref_state=''):
        # Setting the inputs before modal analysis
        self.inputs = inputs
        
        # Running the PowerFactory modal analysis
        # - unnecessarily computing modal analysis just to get the
        #   states; should be fixed in the future
        self.stc.modal_analysis(self._tempdir)

        # Getting the states dictionaries
        states = copy.deepcopy(self.states)
        
        # Getting the dimensions 
        n_x = len(states['states_idx']) # number of states 
        
        # Remove a state in the case it's, e.g., a reference angle
        if ref_state:
            model_name, state_variable = re.match(r'(.+)\\(s:.+)',ref_state).groups()
            states['variables'][model_name].remove(state_variable)
            if not states['variables'][model_name]:
                states['variables'].pop(model_name)
            states['states_idx'].pop(ref_state)
            n_x-=1
        
        # Prepare the outputs used for linearization
        states_outputs = {'ElmRes' : 'ModalComp', 'variables':{}}
        states_outputs['variables'].update(states['variables'])
        states_outputs['variables'].update(outputs['variables'])
        temp_elmres = self.elmres # Used to return to the original outputs after linearization
        self.outputs = states_outputs
        
        # Get a list of ordered states
        ordered_states = [state for state, _ in sorted(states['states_idx'].items(), key = lambda item: item[1])]
        # Get a list of  ordered outputs    
        ordered_outputs = [output for output, _ in sorted(outputs['outputs_idx'].items(), key = lambda item: item[1])]
        
        # Simulate the system with specified inputs
        wave_specs = inputs[0]['wave_specs']
        self.stc.inc.iopt_fastout = 0
        self.stc.inc.dtout = wave_specs['step']
        res = self.simulate(wave_specs['tstart'], wave_specs['tstop'], wave_specs['step'])
        res -= res.iloc[0,:] # Centering with respect to steady-state
       
        #Y = np.vstack((res[ordered_states].iloc[2:].values.T,res[ordered_outputs].iloc[1:-1].values.T/1000))
        #u = np.vstack((wavegen.wavetype[inpt['wave_specs']['type']](**inpt['wave_specs'])['y1'] for inpt in inputs))
        #phi = np.vstack((res[ordered_states].iloc[1:-1].values.T, u[:,:-1]-1))
        
        Y = np.vstack((res[ordered_states].iloc[1:-1].values.T,res[ordered_outputs].iloc[0:-2].values.T))
        inpts = tuple(wavegen.wavetype[inpt['wave_specs']['type']](**inpt['wave_specs'])['y1'] for inpt in inputs)
        u = np.vstack(inpts)
        phi = np.vstack((res[ordered_states].iloc[0:-2].values.T, u[:,:-1]-1))
        
        # Least-squares solution
        theta = np.linalg.solve(phi@phi.T, phi@Y.T).T
        
        # Unpacking theta
        A = theta[:n_x,:n_x]
        B = theta[:n_x,n_x:]
        C = theta[n_x:,:n_x]
        D = theta[n_x:,n_x:]
        self.elmres = temp_elmres
        self.stc.inc.p_resvar = self.elmres.pfobject
        return A, B, C, D

def generate_inputs(inputs, replace=False):
    for inpt in inputs:
        if inpt['source'] == 'generate':
            wave_specs = inpt['wave_specs']
            signal = wavegen.wavetype[wave_specs['type']](**wave_specs)
            filepath = inpt['filepath']
            sigmat = np.empty((signal['time'].shape[0],len(signal))) # Allocate a matrix of size (len(time)) x (number of signals)
            sigmat[:,0] = signal['time'] # Write the time vector
            for signum in range(1,len(signal)): # Iterate over the rest of the signals
                    sigmat[:,signum] = signal['y'+str(signum)] # Write the rest of the signals
            if os.path.isdir(filepath):
                filepath = os.path.join(filepath,
                                        'signal_{}.csv'.format(inpt['ElmFile']))
            with open(filepath, 'w+', newline='') as csvfile: # Write the waveform to the provided filepath
                writer = csv.writer(csvfile, delimiter=' ')
                csvfile.write(str(len(signal)-1)+'\n')
                writer.writerows(sigmat)
            if replace:
                inpt['source'] = 'ext_file'

def d2c(a, b, c ,d, Ts ,method='bi'):
    n=a.shape[0]
    nb=b.shape[1]
    nc=c.shape[0]
    tol=1e-12
    
    if method=='zoh':
            tmp1=np.hstack((a,b))
            tmp2=np.hstack((np.zeros((nb,n)),np.eye(nb)))
            tmp=np.vstack((tmp1,tmp2))
            s=logm(tmp)
            s=s/Ts
            #if norm(np.imag(s),inf) > sqrt(sp.finfo(float).eps):
            #    print("Warning: accuracy may be poor")
            s=np.real(s)
            A=s[0:n,0:n]
            B=s[0:n,n:n+nb]
            C=c
            D=d
    elif method=='foh':
        a=np.mat(a)
        b=np.mat(b)
        c=np.mat(c)
        d=np.mat(d)
        Id = np.mat(np.eye(n))
        A = logm(a)/Ts
        A = np.real(np.around(A,12))
        Amat = np.mat(A)
        B = (a-Id)**(-2)*Amat**2*b*Ts
        B = np.real(np.around(B,12))
        Bmat = np.mat(B)
        C = c
        D = d - C*(Amat**(-2)/Ts*(a-Id)-Amat**(-1))*Bmat
        D = np.real(np.around(D,12))
    elif method=='bi':
        a=np.mat(a)
        b=np.mat(b)
        c=np.mat(c)
        d=np.mat(d)
        poles=eigvals(a)
        #if any(abs(poles-1)<200*sp.finfo(float).eps):
        #    print("d2c: some poles very close to one. May get bad results.")
        
        I=np.mat(np.eye(n,n))
        tk = 2 / np.sqrt (Ts)
        A = (2/Ts)*(a-I)*inv(a+I)
        iab = inv(I+a)*b
        B = tk*iab
        C = tk*(c*inv(I+a))
        D = d- (c*iab)
    else:
        print("Method not supported")
        return
    return A, B, C, D

def start_openmodelica():
    global omchandle
    omchandle = OMCSessionZMQ()
    atexit.register(stop_openmodelica)

def stop_openmodelica():
    global omchandle
    del omchandle