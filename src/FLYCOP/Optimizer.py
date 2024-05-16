# Clase abstracta que define el comportamiento de un optimizador de consortios donde se podrían usar diferentes algoritmos de optimización como los de SMAC3
#, DEAP, etc.
from abc import ABC, abstractmethod

class ParameterOptimizer(ABC):
    @abstractmethod
    def load_configuration_space(self, configuration_space):
        pass

    @abstractmethod
    def run_optimization(self, objective_function):
        pass

    @abstractmethod
    def get_best_parameters(self):
        pass

    @abstractmethod
    def visualize_results(self):
        pass
