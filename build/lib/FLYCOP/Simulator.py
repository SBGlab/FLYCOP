from abc import ABC, abstractmethod

class Simulator(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def load_model(self, model_file):
        pass

    @abstractmethod
    def load_consortia(self,consortia,initial_biomass):
        pass

    @abstractmethod
    def simulate(self, simulation_params):
        pass

    @abstractmethod
    def get_results(self):
        pass

    @abstractmethod
    def visualize_results(self):
        pass
