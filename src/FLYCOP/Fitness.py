# Class to implement a fitness function repository
class Fitness:
    def __init__(self):
        self.fitness_functions = {}
        self.add_function('metabolite_maximization', self.metabolite_maximization)
        self.add_function('biomass_maximization', self.biomass_maximization
        self.default_fitness='metabolite_maximization'

    def add_fitness_function(self, function_name,fitness_function):
        self.fitness_functions[function_name]=fitness_function
        self.default_fitness=function_name

    def remove_fitness_function(self, fitness_function):
        if fitness_function in self.fitness_functions:
            self.fitness_functions.remove(fitness_function)
        else:
            print(f"{fitness_function} is not in the repository.")

    def list_fitness_functions(self):
        print("Fitness functions in the repository:")
        for idx, fitness_function in enumerate(self.fitness_functions, start=1):
            print(f"{idx}. {fitness_function}")

    def set_default_fitness(self,fitness_name):
        if fitness_name in self.fitness_functions:
            self.default_fitness=fitness_name
        else:
            print(f"{fitness_name} is not in the repository.")

    def execute_default_fitness(self,*args,**kwargs):
        return self.fitness_functions[self.default_fitness](*args,**kwargs)

    def execute_fitness_function(self, fitness_name,*args, **kwargs):
        return self.fitness_functions[fitness_name](*args,**kwargs)

    def metabolite_maximization(self,metabolites_conc,metabolite_id):
        # Function to maximize the production of a metabolite
        # Input: Metabolite ID and a pandas dataframe with the metabolite concentrations
        # Output: calculate the metabolite concentration at the last time point minus 1
        # Adds a fitness function to maximize the production of a metabolite

        return(1/metabolites_conc.iloc[-2][metabolite_id])
    def biomass_maximization(self,biomasses):
        # Function to maximize the production of biomass
        # Input: A data.frame with biomass concentrations
        # Output: None
        # Adds a fitness function to maximize the production of biomass
        # Get the biomass column names, the columns names contain the word "biomass"
        columns_biomass = biomasses.filter(like='biomass').columns
        biomasses_values= biomasses.iloc[-2][columnas_biomass]
        return(1/sum(biomasses_values))