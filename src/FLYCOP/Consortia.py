from cobra import Model
from cobra.io import load_model
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

# Class to design a consortium of species
class Consortia:
    name = ""
    species = []

    def __init__(self, name, species=None):
        self.name = name
        if species is None:
            self.species = []
        else:
            # Check if all elements of species are Model objects
            if all(isinstance(specie, Model) for specie in species):
                self.species = species
            else:
                raise TypeError('The species must be a list of cobra.Model objects')

    def uniformize():
        # Function to standardize the models of consortium species with a common namespace for metabolites and reactions
        # Uses the BiGG namespace
        # Input: None
        # Output: None
        pass

    def add_species(self, species_model):
        self.species.append(species_model)

    def add_bigg_species(self, bigg_id):
        # Function to add a species to the consortium based on a BiGG ID
        # Input: BiGG ID of the species
        # Output: None
        # Downloads the model of the species based on its BiGG ID and adds it to the consortium
        model = load_model(bigg_id)
        self.species.append(model)

    def remove_species(self, species_model):
        if species_model in self.species:
            self.species.remove(species_model)
        else:
            print(f"{species_model} is not in the consortia.")

    def list_species(self):
        print(f"Species in {self.name}:")
        for idx, species_model in enumerate(self.species, start=1):
            print(f"{idx}. {species_model}")

    def get_species(self):
        return self.species

    def get_species_names(self):
        return [species.name for species in self.species]

    def get_species_by_name(self, name):
        for species in self.species:
            if species.name == name:
                return species
        return None

    def get_common_metabolites(self, specieName1, specieName2):
        # Function to get common metabolites between two species
        # Input: Names of the two species
        # Output: List of common metabolites
        specie1 = self.get_species_by_name(specieName1)
        specie2 = self.get_species_by_name(specieName2)
        metabolites1 = specie1.metabolites
        metabolites2 = specie2.metabolites
        common_metabolites = []
        for metabolite in metabolites1:
            if metabolite in metabolites2:
                common_metabolites.append(metabolite)
        return common_metabolites

    def get_common_reactions(self, specieName1, specieName2):
        # Function to get common reactions between two species
        # Input: Names of the two species
        # Output: List of common reactions
        specie1 = self.get_species_by_name(specieName1)
        specie2 = self.get_species_by_name(specieName2)
        reactions1 = specie1.reactions
        reactions2 = specie2.reactions
        common_reactions = []
        for reaction in reactions1:
            if reaction in reactions2:
                common_reactions.append(reaction)
        return common_reactions

    def get_consortia_communications(self):
        # Function to obtain communications between consortium species based on exchange metabolites
        # Input: None
        # Output: List of tuples with communications between consortium species
        communications = []
        for idx, specie in enumerate(self.species, start=1):
            for metabolite in specie.metabolites:
                if metabolite.id.endswith("_e") and any(reaction.id.startswith("EX_") for reaction in metabolite.reactions):
                    # Check that the metabolite has an exchange reaction with ID EX_XXXXX
                    for specie2 in self.species[idx:]:
                        if metabolite in specie2.metabolites:
                            communications.append((specie.name, specie2.name, metabolite.id))
        return communications

    def draw_consortia(self):
        # Function to draw the consortium, species are drawn along with communications between them,
        # and exchange metabolites are displayed on the edges representing communications
        # Input: None
        # Output: None
        communications = self.get_consortia_communications()
        G = nx.Graph()
        for specie in self.species:
            G.add_node(specie.name)
        for communication in communications:
            G.add_edge(communication[0], communication[1])
            # Include the exchange metabolite in the communication
            G[communication[0]][communication[1]]['metabolite'] = communication[2]
            pos = graphviz_layout(G, prog="neato")
            nx.draw(G, pos, with_labels=True, arrows=True)
        plt.show()

    def get_consortia_communications_flux(self):
        # Function to obtain communications between consortium species based on exchange metabolites,
        # restricted to those with fluxes compatible with the possible exchange fluxes given by an FVA.
        # Input: None
        # Output: List of tuples with communications between consortium species
        communications = []
        for idx, specie in enumerate(self.species, start=1):
            # Open all exchange reactions of the species by modifying the bounds to (-1000, 1000)
            for reaction in specie.reactions:
                if reaction.id.startswith("EX_"):
                    reaction.bounds = (-1000, 1000)

            specieA_fva = cobra.flux_analysis.flux_variability_analysis(specie, loopless=True)
            # Use FVA to see exchange reactions with fluxes compatible with the possible exchange fluxes of the other species
            # Negative flux in the one that consumes and positive flux in the one that produces
            for metabolite in specie.metabolites:
                if metabolite.id.endswith("_e") and any(reaction.id.startswith("EX_") for reaction in metabolite.reactions):
                    # Check that the metabolite has an exchange reaction with ID EX_XXXXX
                    # Check the flux of the exchange reaction
                    if metabolite in specie2.metabolites:
                        flux_specieA_max = specieA_fva.loc["EX_" + metabolite.id[:-2], "maximum"]
                        flux_specieA_min = specieA_fva.loc["EX_" + metabolite.id[:-2], "minimum"]
                        for specie2 in self.species[idx:]:
                            specieB_fva = cobra.flux_analysis.flux_variability_analysis(specie2, loopless=True)
                            flux_specieB_max = specieA_fva.loc["EX_" + metabolite.id[:-2], "maximum"]
                            flux_specieB_min = specieA_fva.loc["EX_" + metabolite.id[:-2], "minimum"]
                            # Check if the fluxes are compatible, if the maximum or minimum of specieA has the opposite sign to the maximum or minimum of specieB
                            if (flux_specieA_max * flux_specieB_max < 0) or (flux_specieA_min * flux_specieB_min < 0):
                                communications.append((specie.name, specie2.name, metabolite.id))
        return communications

    def get_consortia_design(self, media):
        # Function to obtain communications between consortium species based on exchange metabolites given a specified culture medium
        # For communication search, exchange metabolites of exchange reactions of each species are used along with an FBA simulation
        # to check if there is flux through them
        # Input: Culture medium
        # Output: List of tuples with communications between consortium species
        communications = []
        communications = self.get_consortia_communications(self)

    def draw_consortia(self):
        # Function to draw the consortium, species are drawn along with communications between them,
        # and exchange metabolites are displayed on the edges representing communications
        # Input: None
        # Output: None
        communications = self.get_consortia_communications()
        G = nx.Graph()
        for specie in self.species:
            G.add_node(specie.name)
        for communication in communications:
            G.add_edge(communication[0], communication[1])
            # Include the exchange metabolite in the communication
            G[communication[0]][communication[1]]['metabolite'] = communication[2]
            pos = graphviz_layout(G, prog="neato")
            nx.draw(G, pos, with_labels=True, arrows=True)
        plt.show()


