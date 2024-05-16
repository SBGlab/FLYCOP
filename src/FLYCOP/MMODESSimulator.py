from mmodes import Consortium
from Simulator import Simulator
class MMODESSimulator(Simulator):
    def __init__(self):
        super().__init__()
        self.consortia = Consortium()
        self.initial_abundances=dict()
        self.results = None
        self.simulation_params = None
        self.media=dict()

    def load_model(self, model):
        self.consortia = Consortium(model)

    def load_consortia(self, consortia, initial_biomass):
        for i,model in enumerate(consortia):
            self.consortia.add_model(model, volume_0=initial_biomass[i])
        self.initial_abundances = initial_biomass

    def simulate(self, simulation_params):
        self.simulation_params = simulation_params
        self.consortia.media= self.media
        self.results = self.consortia.run(plot=False)

    def get_results(self):
        return self.results

    def visualize_results(self):
        self.consortia.visualize_results(self.results, self.simulation_params)

    def export_results(self, export_file):
        self.consortia.export_results(self.results, export_file)