from metconsin import surfin_fba, prep_models
from Simulator import Simulator

class SurfinFBASimulator(Simulator):
    def __init__(self):
        super().__init__()
        self.consortia = None
        self.initial_abundances=dict()
        self.results = None
        self.simulation_params = None
        self.media=dict()

    def load_model(self, model):
        self.consortia.add(prep_models.prep_cobrapy_models(model))

    def load_consortia(self, consortia, initial_biomass):
        self.consortia = prep_models.prep_cobrapy_models(consortia)
        self.initial_abundances = initial_biomass

    def simulate(self, simulation_params):
        self.simulation_params = simulation_params
        self.results = surfin_fba(self.consortia, x0=self.initial_abundances,y0=self.media,endtime=self.simulation_params["maxCycles"])

    def get_results(self):
        return self.results

    def visualize_results(self):
        surfin_fba.visualize_results(self.results, self.simulation_params)

    def export_results(self, export_file):
        surfin_fba.export_results(self.results, export_file)