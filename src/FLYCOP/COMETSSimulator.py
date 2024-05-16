# Clase interfaz con COMETS para definir una simulación FBA
from tkinter import N
import cometspy
from Consortia import Consortia
import warnings
from Simulator import Simulator
import pandas as pd
import copy
import matplotlib as plt

class COMETSSimulator(Simulator):
    def __init__(self):
        super().__init__()
        self.consortia = []
        self.simulation_params =  self.initialize_params(None, None)
        self.set_default_parameters()
        self.results = None
        self.simulation = None
        self.layout=cometspy.layout()

    def set_default_parameters(self):
        # Establece los parámetros por defecto de la simulación FBA
        self.simulation_params.set_param('maxCycles',100)

        ## TODO pegar el codigo del FLYCOP 1.5

    def load_model(self, model_file):
        # Carga un modelo SBML
        comets_model=cometspy.model()
        comets_model.load_cobra_model(model_file)
        comets_model.initial_pop = [0, 0, 0.1]
        self.consortia.append()
        self.layout.add_model(comets_model)


    def load_consortia(self,consortia,initial_biomases=None):
        # Carga un consorcio de modelos SBML
        # Transformar la lista de modelos de consortia a modelos COMETS
        self.layout.add_typical_trace_metabolites()
        # Check if initial_biomass is defined
        if initial_biomases is None:
            initial_biomases = dict()
            for model in consortia:
                initial_biomases[model.id] = 0.1

        for model in consortia:
            #chequear si el modelo es de cobra
            if not isinstance(model, cometspy.model):
                comets_model=cometspy.model()
                #self.__ensure_sinks_are_not_exchanges(model) It is not necessary
                comets_model.load_cobra_model(model)
                comets_model.initial_pop = [0, 0, initial_biomases[model.id]]
            else:
                comets_model = model
            self.consortia.append(comets_model)
            self.layout.add_model(comets_model)
            self.layout.update_models()
        

    def load_media(self,media):
        # Chequear que la entrada es un diccionario con los metabolitos y sus concentraciones
        if type(media) is not dict:
            raise TypeError('The media must be a dictionary')
        #Chequear que los metabolitos del medio de cultivo siguen nomenclatura BiGG
        for metabolite in media.keys():
            # chequear que los metabolitos están en almenos uno de modelos de la lista  de consortia
            if not any([metabolite in list(model.get_exchange_metabolites()) for model in self.consortia]):
                warnings.warn('The metabolite '+metabolite+' is not in any of the models of the consortia')
            self.layout.set_specific_metabolite(metabolite,media[metabolite])
        
    def load_bounds(self,model_id, reaction_id,lb,up):
        # Chequear que el modelo existe
        if not any([model_id in model.id for model in self.consortia]):
            raise ValueError('The model '+model_id+' is not in the consortia')
        # Chequear que la reacción existe
        if not any([reaction_id in model.reactions for model in self.consortia]):
            raise ValueError('The reaction '+reaction_id+' is not in the model '+model_id)
        # Chequear que los límites son números
        if not (type(lb) is int or type(lb) is float) or not (type(up) is int or type(up) is float):
            raise TypeError('The lower and upper bounds must be numbers')
        # Cambiar los límites de la reacción del modelo asociado
        #self.layout.set_specific_reaction_bounds(model_id,reaction_id,lb,up)

    def simulate(self, simulation_params=None):
        # Realiza la simulación FBA
        # Check if simulation_params is None to use self.simulation_params
        if simulation_params is not None:
            self.simulation_params = simulation_params
        self.simulation = cometspy.comets(self.layout,self.simulation_params)
        self.simulation.run()
        self.results = self.simulation # TODO: add more results like metabolites concentrations

    def get_results(self):
        # Devuelve los resultados de la simulación FBA
        return self.results

    def visualize_results(self):
        # Visualiza los resultados de la simulación FBA
        self.simulation.plot()

    def change_infite_value_comets(self,cometsfile):
        #read input file
        fin = open(cometsfile, "rt")
        #read file contents to string
        data = fin.read()
        #replace all occurrences of the required string
        data = data.replace('inf', '1000')
        #close the input file
        fin.close()
        #open the input file in write mode
        fin = open(cometsfile, "wt")
        #overrite the input file with the resulting data
        fin.write(data)
        #close the file
        fin.close()

    def __ensure_sinks_are_not_exchanges(self, cobramodel):
        # Check if the sink reactions are exchanges for all of consortia models
        rxn_names = ['sink_PHAg','SK_pqqA_kt_c']
        for rxn in rxn_names:
            cobramodel.reactions.loc[cobramodel.reactions.REACTION_NAMES == rxn, 'EXCH'] = False
            cobramodel.reactions.loc[cobramodel.reactions.REACTION_NAMES == rxn, 'EXCH_IND'] = 0

    def __intersection_exchange_metabolites(self,models):
        # Calculo del medio de varias especies
        all_exchanged_mets = []
        for m in models:
            all_exchanged_mets.extend(m.get_exchange_metabolites())
        all_exchanged_mets = sorted(list(set(list(all_exchanged_mets))))
        return(all_exchanged_mets)

    def __create_media(self,consortia):
        metabolite_list=self.__intersection_exchange_metabolites(consortia)
        media_dict=dict()
        # By default 0 ammount
        for metabolite in metabolite_list:
            media_dict[metabolite]=10.0
        media_dict["ac_e"]=0.0
        media_dict["o2_e"]=1000.0
        media_dict["nh3_e"]=30.0
        media_dict["ca2_e"]=1000
        media_dict["co2_e"]=100
        media_dict["cobalt2_e"]=1000
        media_dict["cu2_e"]=1000
        media_dict["fe2_e"]=100
        media_dict["fe3_e"]=100
        media_dict["h_e"]=100
        media_dict["h2o_e"]=1000
        media_dict["k_e"]=100
        media_dict["mg2_e"]=100
        media_dict["mn2_e"]=100
        media_dict["mobd_e"]=100
        media_dict["na1_e"]=100
        media_dict["no3_e"]=100
        media_dict["photon_e"]=1000
        media_dict["pi_e"]=100
        media_dict["so4_e"]=100
        media_dict["zn2_e"]=1000
        media_dict["btn_e"]=40
        media_dict["cl_e"]=1000
        media_dict["nh4_e"]=40
        media_dict["ni2_e"]=1000

        # Return pandas dataframe
        return(pd.DataFrame(media_dict.items(),columns=['metabolite','init_amount']))
    
    def set_params(self, parameters:cometspy.params):
        self.simulation_params=parameters

    def initialize_params(self,package, globall):
        """Function to initialize the comets' params class
        it can be initialize by submitting two files, one for the package parameters
        and one for the global ones.
        If you don't submit a file, the params class will be initialize with the values stated below
        which have been tested in this simulation"""
    
        if package and globall is not None:
            params = cometspy.params(global_params = globall, package_params= package)
        elif package is not None:
            params = cometspy.params(package_params=package)
        elif globall is not None:
            params = cometspy.params(global_params=globall)
        else:
            params = cometspy.params()
            params.all_params['maxCycles']=1000
            params.all_params['timeStep']=0.1
            params.all_params['spaceWidth']=0.05
            params.all_params['allowCellOverlap']= True
            params.all_params['cellSize'] = 4.3e-13
            params.all_params['deathRate']= 0.00
            params.all_params['dilFactor'] = 1e-2
            params.all_params['dilTime'] = 12
            params.all_params['numRunThreads']= 8
            params.all_params['maxSpaceBiomass']= 100
            params.all_params['defaultVmax']= 20
            params.all_params['defaultKm']= 0.01
            params.all_params['defaultHill']= 1
            params.all_params['showCycleTime']= True
            params.all_params['geneFractionalCost']= 0
            params.all_params['writeTotalBiomassLog']= True
            params.all_params['writeMediaLog']= True
            params.all_params['writeFluxLog']= True
            params.all_params['useLogNameTimeStamp']= False
            params.all_params['flowDiffRate']= 1e-7
            params.all_params['growthDiffRate']=1e-7
            params.all_params['FluxLogRate']= 1
            params.all_params['MediaLogRate']= 1
            params.all_params['numExRxnSubsteps']= 12
            params.all_params['geneFractionalCost'] = 0
            params.all_params['exchangestyle']= 'Standard FBA'
            params.all_params['biomassMotionStyle']= 'Diffusion 2D(Crank-Nicolson)'
            params.all_params['mutRate'] = 1e-9
            params.all_params['addRate'] = 1e-9
            params.all_params['minSpaceBiomass'] = 1e-10

        
        #check if param 'maxSpaceBiomass' has default value
        if (params.all_params['maxSpaceBiomass']== 0.1):
            print('The parameter "maxSpaceBiomass" is set to 0.1.\n' \
                  'It may need to change if the mo used growths well.')

        return params

    def __find_end_cycle(self,simulation_output):
        # it finish when the strains stop the growth
        end_cycle=0
        counter=0
        for index,row in simulation_output.iterrows():
            new_biomasses=row.iloc[1:]
            if index==0 or index==1:
                old_biomasses=row.iloc[1:]
            else:
                result=new_biomasses.subtract(old_biomasses)
                is_growing=result.apply(lambda x: x>0)
                old_biomasses=new_biomasses
                if is_growing.any()==False:
                    counter=counter+1
                    if counter==10:
                        break
            end_cycle=index
        return(end_cycle)

    def __make_df(self, strains, metabolites, comets, max_cycles):
        '''This function creates  a dataframe and the file biomass_vs_met.txt which 
        contais the quantity of each strain and metabolite and has the following columns:
        time(h), strain1 ... strainX, met1 ... metX.'''
        file_name='_'.join(metabolites)
        df = self.simulator.media #We get the media composition results'
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
        #np.savetxt(r'biomass_vs_'+file_name+'_template.txt', df_biomass.values, fmt='%s',delimiter='\t',header='\t'.join(columns)) #The data is saved
        return df_biomass

    def __make_graph(self, df_biomass,strains, metabolites, comets, max_cycles,file_name):
        '''This function creates a figure and saves it to pdf format.
        It also creates the file biomass_vs_met.txt which contais the quantity
        of each strain and metabolite and has the following columns:
        time(h), strain1 ... strainX, met1 ... metX.'''
        file_name='_'.join(metabolites)
        df = comets.media
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
    