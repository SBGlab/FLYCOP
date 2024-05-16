from pickle import NONE
import numpy as np
from ConfigSpace import Configuration, ConfigurationSpace, Float, Integer, Categorical
from matplotlib import pyplot as plt
from Optimizer import ParameterOptimizer
from smac import HyperparameterOptimizationFacade as HPOFacade
from smac import RunHistory, Scenario

class SMAC3Optimizer(ParameterOptimizer):
    def __init__(self, configuration_space=None,objective_function=None,facade=None,parameters={"max_evaluations":100,"deterministic":True}):
        super().__init__()
        self.parameters=parameters
        # Check if the configuration space is None or a ConfigurationSpace class
        if configuration_space is None:
            self.configuration_space=ConfigurationSpace(seed=0)
        if configuration_space is not None:
            if isinstance(configuration_space,ConfigurationSpace):
                self.configuration_space=configuration_space
            elif isinstance(configuration_space,dict):
                self.load_configuration_space(configuration_space)
            else:
                raise ValueError("The selected configuration space is invalid")
        #self.scenario=Scenario(self.configuration_space,deterministic=self.parameters["deterministic"], n_trials=self.parameters["max_evaluations"]) # Create a SMAC scenario TODO add n_trials as parameter
        # check if the objective funcion is a function
        if objective_function is not None:
            if not callable(objective_function):
                raise ValueError("Objective function is not a function")
        else:
            self.objective_function=objective_function
        #Check if the Facade is none or a facade HPOFacade class
        if facade is not None:
            if not isinstance(facade,HPOFacade):
                raise ValueError("The selected Facade is invalid")
        else:
            self.SMAC=facade
        print("The optimizer was created")

    def set_optimization_space(self,cs:ConfigurationSpace):
        # Implementación específica para definir el espacio de optimización en SMAC3
        self.configuration_space=cs

    def set_scenario(self,params):
        self.scenario = Scenario(configspace=self.configuration_space,output_directory="tmp",walltime_limit=params["walltime_limit"],n_trials=params["n_trials"])

    def optimize(self,target):
        self.SMAC = HPOFacade(self.scenario, target,overwrite=True)
        result=self.SMAC.optimize()
        return result


    def load_configuration_space(self, cs):
        # Implementación específica para cargar el espacio de configuración en SMAC3
        conf=ConfigurationSpace(seed=0)
        for x in cs:
            if x['type'] == 'float':
                conf.add_hyperparameters([Float(x['name'], (x['min'], x['max']))])
            elif x['type'] == 'int':
                conf.add_hyperparameters([Integer(x['name'], (x['min'], x['max']))])
            elif x['type'] == 'categorical':
                conf.add_hyperparameters([Categorical(x['name'], x['values'])])
        self.configuration_space=conf

    def load_scenario(self, scenario):
        # Implementación específica para cargar el escenario en SMAC3, la entrada es un diccionario de valores
        #Check if the keys included in the scenario are valid
        
        valid_keys=["n_trials","deterministic"] #TODO add more keys if needed
        for key in scenario.keys():
            if key not in valid_keys:
                raise ValueError(key+"is invalid key in scenario")
        # Check if the configuration space is loaded
        if self.configuration_space is None:
            raise ValueError("Configuration space not loaded")
        else:
            self.scenario=Scenario(configspace=self.configuration_space,n_trials=scenario["n_trials"],deterministic=scenario["deterministic"])

    def load_objective_function(self, objective_function):
        # Implementación específica para cargar la función objetivo en SMAC3
        #Check if the objective function is a function
        if not callable(objective_function):
            raise ValueError("Objective function is not a function")
        self.objective_function=objective_function
        
    def run_optimization(self):
        # Implementación específica para ejecutar la optimización en SMAC3
        # Aquí, 'objective_function' es la función que queremos optimizar
        # Definir el espacio de búsqueda de parámetros
        # Ejecutar la optimización
        #Check if scenario is loaded
        if self.scenario is None:
            raise ValueError("Scenario not loaded")
        #Check if objective function is loaded
        if self.objective_function is None:
            raise ValueError("Objective function not loaded")
        #Check if the configuration space is loaded
        if self.configuration_space is None:
            raise ValueError("Configuration space not loaded")

        # Check if SMAC is loaded
        if self.SMAC is None:
            self.SMAC = HPOFacade(self.scenario, self.objective_function,overwrite=True)

        result=self.SMAC.optimize()
        return(result)



    def get_best_parameters(self):
        # Implementación específica para obtener los mejores parámetros de SMAC3
        pass

    def visualize_results(self):
        # Implementación específica para visualizar los resultados de SMAC3
        # Plot the trajectory of the optimization
        self.SMAC.solver.trajectory.plot_trajectory()
        plt.show()
        # Plot the performance of incumbents
        self.SMAC.solver.trajectory.plot_incumbent_trajectory()
        plt.show()
        # Visualize the performance of all evaluated configurations
        self.SMAC.solver.trajectory.plot_all()
        plt.show()
        # Visualize the optimization trace
        self.SMAC.solver.plot_cost_over_time()
        plt.show()
