import os
import gurobipy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from FLYCOP import FLYCOP
import copy
import cobra
import cometspy as c
from ConfigSpace import Configuration, ConfigurationSpace, Float, Categorical,Integer
from smac import RunHistory, Scenario
import matplotlib
matplotlib.use('Agg') 
def make_df_and_graph(strains, metabolites, comets, max_cycles):
    '''This function creates a figure and saves it to pdf format.
    It also creates the file biomass_vs_met.txt which contais the quantity
    of each strain and metabolite and has the following columns:
    time(h), strain1 ... strainX, met1 ... metX.'''
    file_name='_'.join(metabolites)
    df = comets.media #We get the media composition results'
    df_media=copy.deepcopy(df.loc[df['cycle']<max_cycles])
    df2=comets.total_biomass #We get the biomass results
    df_biomass=copy.deepcopy(df2.loc[df2['cycle']<max_cycles])
    columns=['cycle']
    for i in range(0,len(strains)):
        columns.append(strains[i])
    df_biomass.columns=columns
    """For each metabolite a column with all zeros is added to the dataframe and each row that contains a value
     (metabolite concentration)is changed in the dataframe"""
    for d in metabolites:
        columns.append(d)
        met =df_media.loc[df_media['metabolite'] == d]
        temp=np.zeros(max_cycles) #Create an array with all zeros
        df_biomass[d]=temp #We added it to the dataframe
        j=1
        while j < (max_cycles): #For each cycle
            if (met.cycle==j).any(): #If the row exists
                df_biomass.loc[j-1,d] = float(met[met['cycle']==j]['conc_mmol']) #Its dataframe value is changed
            j+=1
    df_biomass.columns=columns
    np.savetxt(r'biomass_vs_'+file_name+'_template.txt', df_biomass.values, fmt='%s',delimiter='\t',header='\t'.join(columns)) #The data is saved
    #---Starting the figure
    plt.ioff()
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('biomass (g/L)')
    c=['k', 'm', 'b', 'g', 'r']
    j=0
    for i in strains:
        ax.plot(df_biomass['cycle']*0.1, df_biomass[i], label=i, color=c[j])
        j+=1
    ax2 = ax.twinx()
    ax2.set_ylabel('metabolite conc (mM)')
    for m in metabolites:
        ax2.plot(df_biomass['cycle']*0.1, df_biomass[m], label=m)
    handles, labels = ax.get_legend_handles_labels()
    handle_list, label_list = [], []
    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    handles, labels = ax2.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    plt.legend(handle_list, label_list)
    #saving the figure to a pdf
    plt.savefig('biomass_vs_'+file_name+'_template_plot.pdf')
    return df_biomass
###################################################################
#
#               Violacein producer- Model generation
#
###################################################################
#Introduction of needed reactions
def violacein_model():
    import cobra as cb
    model_producer = cb.io.load_model('iWFL_1372')
    # vioA

    reaction = cb.Reaction('vioA')
    reaction.name = 'Flavin-dependent L-tryptophan oxidase'
    reaction.subsystem = ''
    reaction.lower_bound = -1000  #
    reaction.upper_bound = 1000  # Reversible, as directionality is unknown

    #metabolites associated with the reaction

    trp__L_c=model_producer.metabolites.trp__L_c
    o2_c=model_producer.metabolites.o2_c
    h2o2_c=model_producer.metabolites.h2o2_c
    i3i3ppa_c=cb.Metabolite(
        '2i3i3ppa_c',
        formula='C11H9N2O2',
        name='2-imine-3-(indol-3-yl)propanoate',
        compartment='c',
        charge=-1)

    reaction.add_metabolites({
        trp__L_c: -1,
        o2_c: -1,
        h2o2_c: 1,
        i3i3ppa_c:1})

    reaction.gene_reaction_rule = ''

    model_producer.add_reactions([reaction])
    #Introduction of needed reactions
    # vioB
    reaction = cb.Reaction('vioB')
    reaction.name = '2-imino-3-(indol-3yl)propanoate dimerase'
    reaction.subsystem = ''
    reaction.lower_bound = -1000  #
    reaction.upper_bound = 1000  # Reversible, as directionality is unknown

    #metabolites associated with the reaction
    i3i3ppa_c=model_producer.metabolites.get_by_id('2i3i3ppa_c')
    i3pyridm_c=cb.Metabolite(
        'i3pyridm_c',
        formula='C22H18N4O4',
        name='indole-3-pyruvate imine dimer ',
        compartment='c',
        charge=0)

    reaction.add_metabolites({
        i3i3ppa_c: -2,
        i3pyridm_c:1})

    reaction.gene_reaction_rule = ''

    model_producer.add_reactions([reaction])
    #Introduction of needed reactions
    # vioE
    reaction = cb.Reaction('vioE')
    reaction.name = 'Prodeoxyviolacein synthase'
    reaction.subsystem = ''
    reaction.lower_bound = 0  #
    reaction.upper_bound = 1000  # This is the default

    #metabolites associated with the reaction

    i3pyridm_c=model_producer.metabolites.i3pyridm_c
    ptdeovio_c=cb.Metabolite(
        'ptdeovio_c',
        formula='C21H15N3O2',
        name='protodeoxyviolaceinic acid',
        compartment='c',
        charge=0)

    reaction.add_metabolites({
        i3pyridm_c: -1,
        ptdeovio_c:1})

    reaction.gene_reaction_rule = ''

    model_producer.add_reactions([reaction])
    #Introduction of needed reactions
    # vioD
    reaction = cb.Reaction('vioD')
    reaction.name = 'Protodeoxyviolaceinate monooxygenase'
    reaction.subsystem = ''
    reaction.lower_bound = 0  #
    reaction.upper_bound = 1000  # This is the default

    #metabolites associated with the reaction

    ptdeovio_c=model_producer.metabolites.ptdeovio_c
    o2_c=model_producer.metabolites.o2_c
    nadph_c=model_producer.metabolites.nadph_c
    h_c=model_producer.metabolites.h_c
    nadp_c=model_producer.metabolites.nadp_c
    h2o_c=model_producer.metabolites.h2o_c
    ptvio_c=cb.Metabolite(
        'ptvio_c',
        formula='C21H15N3O3',
        name='protoviolaceinic acid',
        compartment='c',
        charge=0)

    reaction.add_metabolites({
        ptdeovio_c: -1,
        o2_c: -1,
        nadph_c: -1,
        h_c:-1,
        nadp_c : 1,
        h2o_c :1,
        ptvio_c: 1})

    reaction.gene_reaction_rule = ''

    model_producer.add_reactions([reaction])
    #Introduction of needed reactions
    # vioC
    reaction = cb.Reaction('vioC')
    reaction.name = 'Violacein synthase'
    reaction.subsystem = ''
    reaction.lower_bound = 0  #
    reaction.upper_bound = 1000  # This is the default

    #metabolites associated with the reaction

    ptvio_c=model_producer.metabolites.ptvio_c
    o2_c=model_producer.metabolites.o2_c
    nadph_c=model_producer.metabolites.nadph_c
    h_c=model_producer.metabolites.h_c
    nadp_c=model_producer.metabolites.nadp_c
    h2o_c=model_producer.metabolites.h2o_c
    violacein_c=cb.Metabolite(
        'violacein_c',
        formula='C20H13N3O3',
        name='Violacein',
        compartment='c',
        charge=0)

    reaction.add_metabolites({
        ptvio_c: -1,
        o2_c: -1,
        nadph_c: -1,
        h_c:-1,
        nadp_c : 1,
        h2o_c :1,
        violacein_c: 1})

    reaction.gene_reaction_rule = ''

    model_producer.add_reactions([reaction])
    model_producer.genes.get_by_id("ECW_m4007")
    model_producer.reactions.SERD_L.bounds=(0,0)
    model_producer.reactions.CYSDS.bounds=(0,0)
    model_producer.reactions.TRPAS2.bounds=(0,0)
    model_producer.reactions.EX_trp__L_e.bounds=(-5,5)
    #Introduction of needed reactions
    # Diffusion of violacein
    reaction = cb.Reaction('VIOtr')
    reaction.name = 'Violacein diffusion'
    reaction.subsystem = ''
    reaction.lower_bound = 0  #
    reaction.upper_bound = 1000  # Violacein can only go outside the cell
    #metabolites associated with the reaction
    violacein_c=model_producer.metabolites.violacein_c
    violacein_e=cb.Metabolite(
        'violacein_e',
        formula='C20H13N3O3',
        name='violacein',
        compartment='e',
        charge=0)
    reaction.add_metabolites({
        violacein_c: -1,
        violacein_e:1})
    reaction.gene_reaction_rule = ''
    model_producer.add_reactions([reaction])
    model_producer.add_boundary(model_producer.metabolites.get_by_id("violacein_e"), type="exchange")
    #model_producer.reactions.vioA.bounds=(0.002694436,1000)
    #model_producer.reactions.EX_gal_e.bounds=(-1.72,0)
        # open the exchanges for the carbon sources
    for metabolite in ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]:
        exchange_name="EX_"+metabolite
        model_producer.reactions.get_by_id(exchange_name).lower_bound=-1.72
    cobra.io.write_sbml_model(model_producer,'violacein_producer.xml')
    return model_producer
###################################################################
#
#               Trp transformer- Model generation
#
###################################################################
def trp_transformer_model():
    model_transformer = cobra.io.load_model('iEC1372_W3110')
    #model_transformer.reactions.TRPAS2.bounds=(-1000,-0.25)
    #model_transformer.reactions.EX_gal_e.bounds=(-1.72,1000)
    model_transformer.reactions.EX_glc__D_e.bounds=(0,1000)
    model_transformer.reactions.EX_trp__L_e.bounds=(0,1000)

    # open the exchanges for the carbon sources
    for metabolite in ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]:
        exchange_name="EX_"+metabolite
        model_transformer.reactions.get_by_id(exchange_name).lower_bound=-1.72
    cobra.io.write_sbml_model(model_transformer,'trp_transformer.xml')
    return model_transformer

def test_violacein_model(model_producer, violacein_flux, carbon_source):
    # vamos a testear el modelo de violacein con unos parametros de consumo de trp y de producción de violacein para ver si funciona correctamente antes de lanzarlo con FLYCOP
    # ponemos todas las fuentes de carbono a 0 menos la que queremos usar en ese momento, para que el modelo no use otras fuentes de carbono que no sean la que queremos analizar en ese momento
    for metabolite in ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]:
        exchange_name="EX_"+metabolite
        if metabolite==carbon_source:
            model_producer.reactions.get_by_id(exchange_name).lower_bound=-1.72
        else:
            model_producer.reactions.get_by_id(exchange_name).lower_bound=0
    # Ahora fijamos el flujo minimo de vioC al valor que queremos testear
    model_producer.reactions.vioC.lower_bound=violacein_flux
    solution=model_producer.optimize()
    print(f"Test of violacein model with carbon source {carbon_source} and violacein flux {violacein_flux}:")
    print(f"Growth rate: {solution.objective_value}, Violacein production: {solution.fluxes['vioC']}")
###################################################################
#
#               COMETS models
#
###################################################################

#converting our models to comets format
model_transformer=trp_transformer_model()
model_producer=violacein_model()
# Vamos a calcular los flujos de trp y de violacein para cada modelo con FVA
# Los calculamos de 10 en 10 entre 0 y 100% y los aplicamos a los bounds de las reacciones correspondientes en cada modelo, para que el espacio de búsqueda de FLYCOP sea más realista. Para ello, se puede crear una función que reciba el modelo, la reacción objetivo (TRPAS2 para el transformador y vioC para el productor) y el porcentaje, y que devuelva los nuevos bounds para esa reacción. Luego, se puede aplicar esa función a cada modelo antes de convertirlo a formato COMETS.
percentages=[0,10,20,30,40,50,60,70,80,90,100]
# Tenemos que guardar los flujos para cada porcentaje y cada fuente de carbono, para luego aplicarlos a los bounds de las reacciones correspondientes en cada modelo antes de convertirlos a formato COMETS.
trp_flux_transformer_dict={}
violacein_flux_producer_dict={}
for carbon_source in ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]:
    for percentage in percentages:
        # tenemos qe poner a 0 los bounds de las fuentes de carbono menos la que estemos usando en ese momento, para que el modelo no use otras fuentes de carbono que no sean la que queremos analizar en ese momento
        for metabolite in ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]:
            exchange_name="EX_"+metabolite
            if metabolite==carbon_source:
                model_transformer.reactions.get_by_id(exchange_name).lower_bound=-1.72
                model_producer.reactions.get_by_id(exchange_name).lower_bound=-1.72
            else:
                model_transformer.reactions.get_by_id(exchange_name).lower_bound=0
                model_producer.reactions.get_by_id(exchange_name).lower_bound=0
        fva_dict_transformer=cobra.flux_analysis.flux_variability_analysis(model_transformer, reaction_list=['TRPAS2'], fraction_of_optimum=percentage/100)
        fva_dict_producer=cobra.flux_analysis.flux_variability_analysis(model_producer, reaction_list=['vioC'], fraction_of_optimum=percentage/100)
        trp_flux_transformer=fva_dict_transformer.loc['TRPAS2','minimum']
        trp_flux_transformer_dict[(carbon_source, percentage)]=trp_flux_transformer
        violacein_flux_producer=fva_dict_producer.loc['vioC','maximum']
        violacein_flux_producer_dict[(carbon_source, percentage)]=violacein_flux_producer
        print(f"Carbon source: {carbon_source}, Percentage: {percentage}%, TRP flux in transformer: {trp_flux_transformer}, Violacein flux in producer: {violacein_flux_producer}")

test_violacein_model(model_producer, violacein_flux_producer_dict[("glc__D_e",80)], "glc__D_e")

model_transformer_comets=c.model(model_transformer)
model_producer_comets=c.model(model_producer)
#saving the models in comets format
model_transformer_comets.id='transformer_strain'
model_producer_comets.id='producer_strain'
###################################################################
#
#               COMETS parameters
#
###################################################################

params=c.params()
params.all_params['maxCycles']=240
params.all_params['timeStep']=0.1
params.all_params['spaceWidth']=0.05
params.all_params['allowCellOverlap']= True
params.all_params['deathRate']= 0.0
params.all_params['maxSpaceBiomass']= 1000
params.all_params['defaultVmax']=20
params.all_params['showCycleTime']=True
params.all_params['useLogNameTimeStamp']=False
params.all_params['FluxLogRate']=1
params.all_params['MediaLogRate']=1
params.all_params['exchangestyle']='Standard FBA'
params.all_params['writeTotalBiomassLog']=True
params.all_params['writeMediaLog']=True

###################################################################
#
#               FLYCOP input space definition
#
###################################################################
# para trp
def violacein_space():
    cs=ConfigurationSpace(seed=0)
    X0=Integer("ratio",bounds=(-1000,1000))
    X1= Categorical("carbon_source",items=["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"])
    X2= Categorical("trp_ratio_flux",items=list([0,10,20,30,40,50,60,70,80,90,100]))
    X3= Categorical("violacein_ratio_flux",items=list([0,10,20,30,40,50,60,70,80,90,100]))
    cs.add_hyperparameters([X0,X1,X2,X3])
    scenario=Scenario(cs, deterministic=True, n_trials=1500)
    return(scenario)
###################################################################
#
#               Scenario definition
#
###################################################################
def carbon_normalization(carbon_source):
    carbon_content={"glc__D_e":6,"gal_e":6,"fru_e":6,"sucr_e":12,"glyc_e":3,"xyl__D_e":5,"mal__D_e":4,"lcts_e":12}
    return carbon_content[carbon_source]

def violacein_scenario(self,parameters,seed: int=0):
    # Creation of the consortium with comets:
    parameters=parameters
    configuration=','.join([str(value) for value in parameters.values()])
    #Leer todos los parametros y sus variables
    #converting our models to comets format
    #layout=c.layout(self.consortia)
    initial_biomass=0.05

    # Trp ratio is the percentage of the trp that goes to the producer, and violacein ratio is the percentage of the trp that goes to violacein production instead of biomass in the producer strain. Both are applied to the initial biomass of the strains, so they affect the initial conditions of the simulation.
    # Con el trp ratio se puede controlar la cantidad de trp que se produce en el transformador y que se consume en el productor, y con el violacein ratio se puede controlar la cantidad de trp que se destina a la producción de violaceina en el productor, lo que afecta directamente a la producción de violaceina. Ambos parámetros son importantes para optimizar la producción de violaceina, ya que un exceso de trp en el transformador puede ser tóxico para las células, y un exceso de trp destinado a la producción de violaceina puede reducir el crecimiento del productor y por tanto su capacidad de producir violaceina a largo plazo.
    # Actualizamos los bounds de las reacciones correspondientes a cada modelo con los valores calculados en el FVA para cada porcentaje, para que el espacio de búsqueda de FLYCOP sea más realista.
    self.simulator.layout.models[0].change_bounds("TRPAS2" , -1000,trp_flux_transformer_dict[(parameters["carbon_source"],parameters["trp_ratio_flux"])])
    self.simulator.layout.models[1].change_bounds("vioC" , 0,violacein_flux_producer_dict[(parameters["carbon_source"],parameters["violacein_ratio_flux"])])
    self.simulator.layout.update_models()

    if parameters["ratio"]>=0:
        biomass_transformer=initial_biomass*parameters["ratio"]
        biomass_producer=initial_biomass
    else:
        biomass_transformer=initial_biomass
        biomass_producer=initial_biomass*(-1*parameters["ratio"])
        
    self.simulator.layout.initial_pop=[[0.0,0.0,biomass_transformer,biomass_producer]] #Initial biomases follow the 99:1 ratio from the paper
    self.simulator.layout.add_typical_trace_metabolites()
    self.simulator.layout.set_specific_metabolite(parameters["carbon_source"],27.7*(6/carbon_normalization(parameters["carbon_source"]))) # normalizado a glucosa
    self.simulator.simulate()
    simulation=self.simulator.get_results()
    print(simulation)
    results=make_df_and_graph([model_transformer_comets.id,model_producer_comets.id],["violacein_e"],simulation,170)
    violacein=results.loc[165,"violacein_e"]
    biomass1=results.iloc[165,1]
    biomass2=results.iloc[165,2]
    if violacein==0:
        fitness=10000
    else:
        fitness=1/violacein
    # Print the parameters and the violacein production
    print("Parameters: "+configuration+",ET: "+str(biomass1)+",EV "+str(biomass2)+", --> Violacein: "+str(violacein)+",Fitness: "+str(fitness))
    print("The triptophane flux in the transformer is: "+str(trp_flux_transformer_dict[(parameters["carbon_source"],parameters["trp_ratio_flux"])]))
    print("The violacein flux in the producer is: "+str(violacein_flux_producer_dict[(parameters["carbon_source"],parameters["violacein_ratio_flux"])]))

    return fitness,0

###################################################################
#
#               FLYCOP run
#
###################################################################

# Consortia list
consortia=[model_transformer_comets,model_producer_comets]
# FLYCOP2 object creation
FLYCOP_optimization=FLYCOP()
# FLYCOP2 set the dynamic FBA selected tool
FLYCOP_optimization.create_simulator(simulator_type="COMETS")
# FLYCOP2 set the input models
FLYCOP_optimization.load_consortia(consortia)
# FLYCOP2 set the COMETS parameters
FLYCOP_optimization.set_simulator_params(params)
# Set the COMETS scenario
FLYCOP_optimization.create_scenario()
FLYCOP_optimization.set_scenario(violacein_scenario)
# Vamos a testear el escenario con unos parámetros concretos para ver si funciona correctamente antes de lanzarlo con SMAC3
test_parameters={"ratio":0.5,"carbon_source":"glc__D_e","trp_ratio_flux":80,"violacein_ratio_flux":80}
fitness, _ = violacein_scenario(FLYCOP_optimization, test_parameters)
# Set the SMAC3 space and run
espace=violacein_space()
result=FLYCOP_optimization.optimize(espace)
# salvamos el resultado
import pickle
with open('violacein_optimization_result.pkl', 'wb') as f:
    pickle.dump(result, f)
    print("Resultado guardado en 'violacein_optimization_result.pkl'")
