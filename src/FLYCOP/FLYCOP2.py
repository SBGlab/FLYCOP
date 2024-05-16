# Main of FLYCOP2

# FLYCOP2 class with the main functions, it contains a consortia object, a set of parameters to optimize, an objective function, 
# a parameter optimizator obj and a dynamic FBA simulation obj called simulator
from tkinter import SE
from SimulatorFactory import SimulatorFactory
from OptimizerFactory import OptimizerFactory
from ConfigSpace import Configuration, ConfigurationSpace, Float
#import Fitness
class FLYCOP2:
    # Constructor
    def __init__(self,simulator_type="COMETS",optimizator_type="SMAC3"):
        self.simulator = self.create_simulator(simulator_type)
        self.parameters = None
        #self.objective_function = Fitness()
        self.scenario=FLYCOP_scenario(self.simulator)
        self.optimizator = self.create_optimizator(optimizator_type)

    def load_consortia(self,model_list=list):
        # Load consortia from mode list
        self.simulator.load_consortia(model_list)

    def create_simulator(self,simulator_type=str):
       # Create a simulator object
        simFactory=SimulatorFactory()
        self.simulator=simFactory.createSimulator(simulator_type)
    
    def create_scenario(self):
        # Create a scenario object
        # Check if the simulator is set
        if self.simulator is None:
            print("Simulator not set")
        else:
            self.scenario=FLYCOP_scenario(self.simulator)

    def set_simulator_params(self,parameters):
        # Set the parameters for the simulator
        self.simulator.set_params(parameters)

    def set_simulator(simulator):
        # Set the simulator object
        self.simulator=simulator

    def simulate(self,parameters):
        # Run the simulation
        self.simulator.simulate(parameters)

    def create_optimizator(self,optimizator_type):
       # Create an optimizer object
        self.optimizator=OptimizerFactory.createOptimizer(optimizer_name=optimizator_type)

    def set_fitness(self,function_name,fitness_function):
       # Set the fitness function
        self.objective_function.add_fitness_function(function_name,fitness_function)

    def set_scenario(self,scenario):
        # Set the scenario, it is an object of the class FLYCOP_scenario
        self.scenario.set_scenario(scenario)
       
    def load_parameters(self,parameters=dict):
        # Load a dictionary of parameters
        self.parameters = parameters

    def optimize(self,scenario):
        self.create_optimizator(optimizator_type="SMAC3")
        #self.optimizator.set_optimization_space(configspace)
        self.optimizator.scenario=scenario
        #self.optimizator.set_scenario({"walltime_limit":120,"n_trials":5})
        result=self.optimizator.optimize(self.scenario.scenario)
        return(result)

    def run_optimization(self):
    # Run the optimization
    # Input: None
    # Output: None
    # Runs the optimization process 
    #Check if the simulator and the optimizator are set
        if self.simulator is None:
                print("Simulator not set")
                return
        if self.optimizator is None:
            print("Optimizator not set")
            return
        #Check if the parameters are set
        if self.parameters is None:
            print("Parameters not set")
            return
        #Check if the objective function is set
        if self.objective_function is None:
            print("Objective function not set")
            return
    #Run the optimization
        self.optimizator.load_objective_function(self.scenario.optimize)
        self.optimizator.run_optimization(self.simulator,self.parameters,self.objective_function)

class FLYCOP_scenario:
# Class to implement a FLYCOP scenario
    def __init__(self,simulator,fitness=None):
        self.simulator=simulator
        self.fitness=fitness
    def configspace(self,configuration_dict) -> ConfigurationSpace:
        # Function to define the configuration space for the optimization
        # Input: None
        # Output: Configuration space
        cs=ConfigurationSpace(configuration_dict)
        return cs

    def scenario(self,parameters):
        pass

    def set_scenario(self,scenario):
        # Function to set the scenario
        # Input: scenario
        # Output: None
        # Asign the new scenario to the object
        setattr(self, "scenario", scenario.__get__(self, FLYCOP_scenario))
        print("The new scenario was created")

    def optimize(self,config:Configuration)->float:
        # Function to run a FLYCOP scenario
        # Input: Configuration
        # Output: fitness value
        
        # Retrieve the parameters to be optimized, biomasses, uptakes, etc.
        for key, value in config.items():
            print(f"{key}: {value}")
            locals()[key] = value
        # Run the simulation
        biomass_variables = [value for key, value in locals().items() if key.startswith('biomass')]
        # Set the biomass variables if they are included in the configuration
        if biomass_variables:
            self.simulator.load_consortia(self.consortia,biomass_variables)
        else:
            self.simulator.load_consortia(self.consortia)
        # Set the bounds parameters if they are included in the configuration they follow the format XXXX_lb and XXXX_ub where XXXX is the name of the reaction
        bounds_variables = {key: value for key, value in locals().items() if key.endswith('_lb') or key.endswith('_ub')}
        # Check if bounds_variables is empty
        if bounds_variables:
            #Extract the reaction name from the key
            bounds_reaction_name= [key.split('_')[0] for key, value in bounds_variables.items()]
            self.simulator.load_bounds(bounds_variables) #TODO create a function to load bounds

        # Set the media variables if they are included in the configuration
        # The media variables follow the format XXXX_e where XXXX is the metabolite name and they dont start with EX_
        media_variables = {key: value for key, value in locals().items() if not key.startswith('EX_') and key.endswith('_e')}
        # Check if media_variables is empty
        if media_variables:
            self.simulator.load_media(media_variables)
        #Run the simulation
        self.simulator.simulate()
        result=self.simulator.get_results()
        #Calculate the fitness value
        fitness_value=self.fitness.execute_fitness_function()
        # Retrieve the fitness value
        return fitness_value


# Function to launch FLYCOP GUI
def launch_GUI():
    pass


# Function to launch FLYCOP command line
def launch_CLI():
    pass

