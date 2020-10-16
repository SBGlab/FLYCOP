#!/usr/bin/python3

############ FLYCOP ############
# Author: Beatriz García-Jiménez
# April 2018
################################

# Example: >>%run ecoliLongTermFLYCOP
#          >>avgfitness,sdfitness=ecoliLongTermFLYCOP_oneConf(-10,-10,-20,-5,-15,-10)
# Goal: individual test to improve consortium {E.coli-E.coli}, depending on glucose, acetate and oxygen uptakes in each of the two strains:
# Run through the function ecoliLongTermFLYCOP_oneConf

 
import cobra #as cb
import pandas as pd
import tabulate
import re
import sys
import getopt
import os.path
import copy
import csv
import math
import cobra.flux_analysis.variability
import massedit
import subprocess
import shutil, errno
import statistics
from cobra import Reaction
#import comets as c


############## Esta función no haría falta ##################################################
### FUNCTION initialize_models #################################    
"""def initialize_models():
 # Only to run 1st time, to build the models!!
 if not(os.path.exists('ModelsInput/iJO1366.mat')):
     print('ERROR! Not iJO1366.mat files with GEM of consortium strains in ModelsInput!')
 else:
  path=os.getcwd()
  os.chdir('ModelsInput')
  model=cobra.io.load_matlab_model('iJO1366.mat')
  # To change the core to the WT biomass as objective
  model.objective="BIOMASS_Ec_iJO1366_WT_53p95M"
  # To remove one copy from GLCtex
  model.reactions.get_by_id('GLCtex_copy2').bounds=(0,0)
  model.reactions.get_by_id('PFK_3').bounds=(0,0) # To avoid to deviate flux of glycolysis, going it through PFK rxn.
  # Replace brackets with compartment location (e.g. "[c]") in metabolite ids by '_' (e.g. "_c") 
  for metabolite in model.metabolites:
    metabolite.id = re.sub('_c$',r'[c]',metabolite.id)
    metabolite.id = re.sub('_p$',r'[p]',metabolite.id)
    metabolite.id = re.sub('_e$',r'[e]',metabolite.id)
    metabolite.id = re.sub('__',r'_',metabolite.id)
    metabolite.compartment = ''
  # To solve possible problems in changing names     
  model.repair()
  # Replace brackets with compartment location (e.g. "[c]") in metabolite ids by '_' (e.g. "_c") 
  for rxn in model.reactions:
    rxn.id = re.sub('_p$',r'(p)',rxn.id)
    rxn.id = re.sub('_c$',r'(c)',rxn.id)
    rxn.id = re.sub('_e$',r'(e)',rxn.id)    
  # To solve possible problems in changing names     
  model.repair()
  cobra.io.save_matlab_model(model,"iJO1366py.mat",'model')
  #
  model.reactions.get_by_id('EX_o2(e)').bounds=(-14.4,1000) # To adjust WT model to realistic values of acetate generation, according to [Steinsiek,2012; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3486086/], with: glc~=10, ac~=3.7, GR~=0.8
  model.reactions.get_by_id('EX_ac(e)').bounds=(-1000,60) # Open lower limit of acetate (allowing to taking it (for cases where there isn't glucose), and limiting the secretion to 3*maximum glucose in our long term experiment (20))
  cobra.io.save_matlab_model(model,"iJO1366py_tmp.mat",'model')
  del(model)
  os.chdir(path)
# end-def"""
################################################################
#loading a sbml model
def initialize_models():
    if not(os.path.exists('ModelsInput/NOMBREDELARCHIVO.sbml')):
        print('ERROR! Not NOMBREDELARCHIVO.sbml files with GEM of consortium strains in ModelsInput!')
    else:
    working_dir=os.getcwd()
    os.chdir('ModelsInput')
    cobramodel = cb.io.read_sbml_model('NOMBREDELARCHIVO.sbml')
    model=c.model(cobramodel)
    #To change the core to the WT biomass as objective
    model.objective="BIOMASS_Ec_iJO1366_WT_53p95M"
    # To remove one copy from GLCtex
    model.change_bounds('GLCtex_copy2',0,0)
    model.change_bounds('PFK_3',0,0) # To avoid to deviate flux of glycolysis, going it through PFK rxn.
    #No need to add an empty location, in layout class we'd have local_media with just a location
    model.write_comets_model(working_dir=working_dir) #Not sure if necessary
    #
    model.change_bounds('EX_o2_e',-14.4,1000) # To adjust WT model to realistic values of acetate generation, according to [Steinsiek,2012; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3486086/], with: glc~=10, ac~=3.7, GR~=0.8
    model.change_bounds('EX_ac_e',-1000,60) # Open lower limit of acetate (allowing to taking it (for cases where there isn't glucose), and limiting the secretion to 3*maximum glucose in our long term experiment (20))
    model.id = model.id + '_tmp'
    model.write_comets_model(working_dir=working_dir) #Not sure if necessary
    model_id=model.id
    return model, model_id #not sure if removing model and read it in initialize_layout, just like in ecolilongtermFlYCOP...

#Stablish the layout file 
def initialize_layout(input_obj):
    layout=c.layout(input_obj)
    #assume it is a layout file, need to add the model with add_models
    #maybe read it with read_comets_model
    layout.add_model(model) #this adds the models to the layout and update the layout with its information, as initial_population
    return layout
    
"""################################################################    
### FUNCTION mat_to_comets #####################################    
# mat_to_comets(modelPath)
# Re-code in python from COMETS code
def mat_to_comets(matInputFile):
    model=cobra.io.load_matlab_model(matInputFile)
    # Open output file:
    with open(matInputFile+'.cmt', mode='w') as f:
        # Print the S matrix
        f.write("SMATRIX  "+str(len(model.metabolites))+"  "+str(len(model.reactions))+"\n")
        for x in range(len(model.metabolites)):
            for y in range(len(model.reactions)):
                if (model.metabolites[x] in model.reactions[y].metabolites):
                    coeff=model.reactions[y].get_coefficient(model.metabolites[x])
                    f.write("    "+str(x+1)+"   "+str(y+1)+"   "+str(coeff)+"\n")
        f.write("//\n")
        
        # Print the bounds
        f.write("BOUNDS  -1000  1000\n");
        for y in range(len(model.reactions)):
            lb=model.reactions[y].lower_bound
            up=model.reactions[y].upper_bound
            f.write("    "+str(y+1)+"   "+str(lb)+"   "+str(up)+"\n")
        f.write("//\n")
        
        # Print the objective reaction
        f.write('OBJECTIVE\n')
        for y in range(len(model.reactions)):
            if (model.reactions[y] in model.objective):
                indexObj=y+1
        f.write("    "+str(indexObj)+"\n")
        f.write("//\n")
        
        # Print metabolite names
        f.write("METABOLITE_NAMES\n")
        for x in range(len(model.metabolites)):
            f.write("    "+model.metabolites[x].id+"\n")
        f.write("//\n")

        # Print reaction names
        f.write("REACTION_NAMES\n")
        for y in range(len(model.reactions)):
            f.write("    "+model.reactions[y].id+"\n")
        f.write("//\n")

        # Print exchange reactions
        f.write("EXCHANGE_REACTIONS\n")
        for y in range(len(model.reactions)):
            if (model.reactions[y].id.find('EX_')==0):
                f.write(" "+str(y+1))
        f.write("\n//\n")                
### end-function-mat_to_comets    
################################################################
"""
def ecolilongTermFLYCOP_oneConf(glu1,ac1,o21,glu2,ac2,o22, fitFunc='MaxYield_MinTime',dirPlot='',repeat=10):
    
    #check if class model exist
   
    try:
        model = model()
    except NameError:
        #clas model doesn't exists
        initialize_models()

    # Determine initial biomasses.
    biomass1=0.01
    biomass2=0.01
    print("Fitness function:"+fitFunc)
    
  # Single GEMs parameter modifications
  # ===================================
    if not(os.path.exists('ecoli_1_tmp.cmd')):
    # 1.1.- [COBRApy] Establish modifications in model 1 
    model1=c.model.read_comets_model('ModelsInput/' + model_id + '.cmd') #loads model from initialize_model
    model1.change_bounds('EX_glc__D_e',-1000,1000)
    model1.change_bounds('EX_o2_e',-1000,1000)
    model1.change_bounds('GLCtex_copy1',0,-glu1)
    model1.change_bounds('O2tex',0,-o21)
    if(ac1<=0):
        model1.change_bounds('ACtex',-1000,-ac1)
        model1.change_bounds('EX_ac_e',ac1,1000)
    else:
        model1.change_bounds('ACtex',-ac1,-1000)
        model1.change_bounds('EX_ac_e',-1000,ac1)
    model1.id='ecoli_1_tmp'
    model1.initial_population=biomass1
    model1.write_comets_model() #save it in the working directory ?
  
    # 1.2.- [COBRApy] Establish modifications in model 2
    model2=c.model.read_comets_model('ModelsInput/' + model_id + '.cmd') #loads model from initialize_model
    model2.change_bounds('EX_glc__D_e',-1000,1000)
    model2.change_bounds('EX_o2_e',-1000,1000)
    model2.change_bounds('GLCtex_copy1',0,-glu2)
    model2.change_bounds('O2tex',0,-o22)
    if(ac1<=0):
        model2.change_bounds('ACtex',-1000,-ac2)
        model2.change_bounds('EX_ac_e',ac2,1000)
    else:
        model2.change_bounds('ACtex',-ac2,-1000)
        model2.change_bounds('EX_ac_e',-1000,ac2)
    model2.id='ecoli_2_tmp'
    model2.initial_population=biomass2
    model2.write_comets_model() #save it in the working directory ?
    
################################################################
### FUNCTION ecoliLongTermFLYCOP_oneConf #######################   
def ecoliLongTermFLYCOP_oneConf(glu1,ac1,o21,glu2,ac2,o22,fitFunc='MaxYield_MinTime',dirPlot='',repeat=10):
  '''
  Call: avgFitness, sdFitness = ecoliLongTerm_oneConf(glu1,ac1,o21,glu2,ac2,o22)

  INPUTS: glu1: lower bound of glucose uptake in model 1.
          ace1: lower bound of acetate uptake in model 1.
          o21: lower bound of oxygen uptake in model 1.
          glu2: lower bound of glucose uptake in model 2.
          ace2: lower bound of acetate uptake in model 2.
          o22: lower bound of oxygen uptake in model 2.
          dirPlot: copy of the graphs with several run results.
          repeat: number of runs with the same configuration.
  OUTPUT: avgFitness: average fitness of 'repeat' COMETS runs with the same configuration (due to it is not deterministic)
          sdFitness: standard deviation of fitness during 'repeat' COMETS runs (see above)
  '''

  if not(os.path.exists('ModelsInput/iJO1366py_tmp.mat')):
      initialize_models()

  # Determine initial biomasses.
  biomass1=0.01
  biomass2=0.01

  print("Fitness function:"+fitFunc)

  # Single GEMs parameter modifications
  # ===================================
    if not(os.path.exists('ecoli_1_tmp.mat.cmt')):
    # 1.1.- [COBRApy] Establish modifications in model 1 
    model1=cobra.io.load_matlab_model('ModelsInput/iJO1366py_tmp.mat')
    model1.reactions.get_by_id('EX_glc__D(e)').bounds=(-1000,1000)
    model1.reactions.get_by_id('EX_o2(e)').bounds=(-1000,1000)
    model1.reactions.get_by_id('GLCtex_copy1').bounds=(0,-glu1)
    model1.reactions.get_by_id('O2tex').bounds=(0,-o21)
    if(ac1<=0):
        model1.reactions.get_by_id('ACtex').lower_bound=-1000
        model1.reactions.get_by_id('ACtex').upper_bound=-ac1
        model1.reactions.get_by_id('EX_ac(e)').bounds=(ac1,1000)
    else:
        model1.reactions.get_by_id('ACtex').lower_bound=-ac1
        model1.reactions.get_by_id('ACtex').upper_bound=1000
        model1.reactions.get_by_id('EX_ac(e)').bounds=(-1000,ac1)
    model1.initial_pop=biomass1
    cobra.io.save_matlab_model(model1,'ecoli_1_tmp.mat','model')
    
      
    # 1.2.- [COBRApy] Establish modifications in model 2
    model2=cobra.io.load_matlab_model('ModelsInput/iJO1366py_tmp.mat')
    model2.reactions.get_by_id('EX_glc__D(e)').bounds=(-1000,1000)
    model2.reactions.get_by_id('EX_o2(e)').bounds=(-1000,1000)
    model2.reactions.get_by_id('GLCtex_copy1').bounds=(0,-glu2)
    model2.reactions.get_by_id('O2tex').bounds=(0,-o22)
    if(ac2<=0):
        model2.reactions.get_by_id('ACtex').lower_bound=-1000
        model2.reactions.get_by_id('ACtex').upper_bound=-ac2
        model2.reactions.get_by_id('EX_ac(e)').bounds=(ac2,1000)
    else:
        model2.reactions.get_by_id('ACtex').lower_bound=-ac2
        model2.reactions.get_by_id('ACtex').upper_bound=1000
        model2.reactions.get_by_id('EX_ac(e)').bounds=(-1000,ac2)
        model2.initial_pop=biomass2

    cobra.io.save_matlab_model(model2,'ecoli_2_tmp.mat','model')
    del(model)
    
    # 2.- [python] 
    mat_to_comets('ecoli_1_tmp.mat')
    mat_to_comets('ecoli_2_tmp.mat')

    # Community parameter modifications
    # =================================        
    # 4.- [shell script] Write automatically the COMETS parameter about initial biomass of strains.
    massedit.edit_files(['ecoliLongTerm_layout_template.txt'],["re.sub(r'XXX','"+str(biomass1)+"',line)"], dry_run=False)
    massedit.edit_files(['ecoliLongTerm_layout_template.txt'],["re.sub(r'YYY','"+str(biomass2)+"',line)"], dry_run=False)        

  # 5.- [COMETS by command line] Run COMETS
  if not(os.path.exists('IndividualRunsResults')):
    os.makedirs('IndividualRunsResults')
  totfitness=0
  sumTotBiomass=0
  sumTotYield=0
  fitnessList=[]
  # To repeat X times, due to random behaviour in COMETS:
  for i in range(repeat):
        with open("output.txt", "w") as f:
            subprocess.call(['./comets_scr','comets_script_template'], stdout=f)
    
        # 6.- [R call] Run script to generate one graph: strains versus metabolite/s
        subprocess.call(['../../Scripts/plot_biomassX2_vs_2mediaItem.sh','template','glc_D','ac','Ecoli1','Ecoli2'])

        # 7.- Compute fitness (measure to optimize):
        print('computing fitness...')
        # 7.1.- Determine endCycle: when glucose and acetate are exhausted
        with open("biomass_vs_glc_D_ac_template.txt", "r") as sources:
            lines = sources.readlines()
            iniPointV=lines[0].split()
            iniBiomass=float(iniPointV[1])+float(iniPointV[2])
            totGlc=float(iniPointV[3])
            endGlcCycle=0
            for line in lines:
                endCycle=int(line.split()[0])
                glcConc=float(line.split()[3])
                acConc=float(line.split()[4])
                if((endGlcCycle==0)and(glcConc==0.0)):
                    endGlcCycle=endCycle
                if((glcConc==0.0)and(acConc==0.0)):
                    break;
                if((glcConc==0.0)and(ac1>=0)and(ac2>=0)):
                    break;
            endPointV=lines[endCycle].split()
        # 7.2.- Compute first element fitness: maximize biomass yield
        # biomass yield= sum(increment in biomass per strain (i.e. biomass final point-biomass initial point))/initial concentration of glucose in the media (total glucose, because the end of our experiment is after glucose finished). In gDW/mmol.
        # To compute final biomass as the maximum biomass of each strain
        finalBiomass1=0
        finalBiomass2=0
        count=0
        for line in lines:
            if(float(line.split()[1])>finalBiomass1):
                finalBiomass1=float(line.split()[1])
            if(float(line.split()[2])>finalBiomass2):
                finalBiomass2=float(line.split()[2])
            if(count>endCycle):
                break;
            count=count+1
        finalBiomass=finalBiomass1+finalBiomass2
        biomassYieldNew=float((finalBiomass-iniBiomass)/(totGlc*0.1801559)) # molecular weigth glucose per mmol
        # For normalizing yield
        MaximumYield=0.6
        # 7.3.- Compute second element fitnes: minimize time        
        fitTime=1-(float(endCycle)/float(240))
        # 7.4.- Compute joint fitness, as a 50% each element.
        if(fitFunc=='Yield'):
            fitness=(biomassYieldNew/MaximumYield) # Normalizing Yield
        elif(fitFunc=='MaxYield_MinTime'):
            fitness=0.5*(biomassYieldNew/MaximumYield)+0.5*fitTime #Normalizing yield
        elif(fitFunc=='YieldNewScattered'): # (biomass^4)*10: To spread values from ~0.45-0.55 values to 0.5 to 1
            fitness=(biomassYieldNew**4)*10
        elif(fitFunc=='MaxYieldNewScattered_MinTime'):
            fitness=0.5*((biomassYieldNew**4)*10)+0.5*fitTime
        elif(fitFunc=='Biomass'):
            fitness=float(finalBiomass-iniBiomass)
        elif(fitFunc=='MaxBiomass_MinTime'):
            fitness=0.5*(float(finalBiomass-iniBiomass))+0.5*fitTime
        elif((fitFunc=='GR')or(fitFunc=='MaxGR_MinTime')):
            numRxnGR1=int(subprocess.check_output(['egrep -A1 "OBJECTIVE" ecoli_1_tmp.mat.cmt | tail -1 | tr -d \[:space:\]'], shell=True))
            numRxnGR2=int(subprocess.check_output(['egrep -A1 "OBJECTIVE" ecoli_2_tmp.mat.cmt | tail -1 | tr -d \[:space:\]'], shell=True))            
            try:
                GR1=float(subprocess.check_output(['egrep "fluxes\{"'+str(endGlcCycle-1)+'"\}\{1\}\{1\}\{1\}" flux_log_template.txt | cut -d"=" -f2 | cut -d" " -f'+str(numRxnGR1+1)], shell=True))
            except:
                GR1=0.0
            try:
                GR2=float(subprocess.check_output(['egrep "fluxes\{"'+str(endGlcCycle-1)+'"\}\{1\}\{1\}\{2\}" flux_log_template.txt | cut -d"=" -f2 | cut -d" " -f'+str(numRxnGR2+1)], shell=True))
            except:
                GR2=0.0
            fitGR=(GR1+GR2)/2
            if(fitFunc=='GR'):
                fitness=fitGR
            elif(fitFunc=='MaxGR_MinTime'):
                fitness=0.5*fitGR+0.5*fitTime

        # To avoid unrealistic cases, because with 10mM of glc the strains can't reach more than ~1 gr/L. I'm not sure if the relation is lineal with less or more glucose, but this solution is better than >1, which it will be very ad-hoc to totGlc=10.
        #if(float(finalBiomass-iniBiomass) > (totGlc/10)):
        if(float(finalBiomass-iniBiomass) > 1.03):  # Given that with both strains with WT, total biomass=1.028
            fitness=0

        # Compute acetate uptake
        numRxnExAc=37 # Position EX_ac(e) in .mat.cmt - Position first rxn + 1 
        # flux in cycle 2 (the first one is usually 0)
        uptakeAc1=float(subprocess.check_output(['egrep "fluxes\{2\}\{1\}\{1\}\{1\}" flux_log_template.txt | cut -d"=" -f2 | cut -d" " -f'+str(numRxnExAc)], shell=True))
        uptakeAc2=float(subprocess.check_output(['egrep "fluxes\{2\}\{1\}\{1\}\{2\}" flux_log_template.txt | cut -d"=" -f2 | cut -d" " -f'+str(numRxnExAc)], shell=True))
            
        print(" Total biomass: "+str(round(finalBiomass,6))+" in cycle "+str(endCycle)+". Biomass yield="+str(round(biomassYieldNew,6)))

        totfitness=totfitness+fitness
        fitnessList.append(fitness)
        sumTotBiomass=sumTotBiomass+finalBiomass
        sumTotYield=sumTotYield+biomassYieldNew
        
        # Copy individual solution
        file='IndividualRunsResults/'+'biomass_vs_glc_D_ac_run'+str(i)+'_'+str(fitness)+'_'+str(endCycle)+'.pdf'        
        shutil.move('biomass_vs_glc_D_ac_template_plot.pdf',file)
        if(dirPlot != ''):
            file2=dirPlot+'biomass_vs_glc_D_ac_'+str(glu1)+'_'+str(ac1)+'_'+str(o21)+'_'+str(glu2)+'_'+str(ac2)+'_'+str(o22)+'_'+str(round(uptakeAc1,1))+'_'+str(round(uptakeAc2,1))+'_run'+str(i)+'_'+str(fitness)+'_'+str(endCycle)+'.pdf'
            shutil.copy(file,file2)
        file='IndividualRunsResults/'+'total_biomass_log_run'+str(i)+'.txt'
        shutil.move('total_biomass_log_template.txt',file)
        file='IndividualRunsResults/'+'media_log_run'+str(i)+'.txt'
        shutil.move('media_log_template.txt',file)
        file='IndividualRunsResults/'+'flux_log_run'+str(i)+'.txt'
        shutil.move('flux_log_template.txt',file)   
    
  avgfitness=totfitness/repeat
  sdfitness=statistics.stdev(fitnessList)
  avgBiomass=sumTotBiomass/repeat
  avgYield=sumTotYield/repeat
  print("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle")
  print(fitFunc+"\t"+str(glu1)+','+str(ac1)+','+str(o21)+','+str(glu2)+','+str(ac2)+','+str(o22)+','+str(round(uptakeAc1,1))+','+str(round(uptakeAc2,1))+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle))
  with open(dirPlot+"configurationsResults"+fitFunc+".txt", "a") as myfile:
      myfile.write("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle\n")
      myfile.write(fitFunc+"\t"+str(glu1)+','+str(ac1)+','+str(o21)+','+str(glu2)+','+str(ac2)+','+str(o22)+','+str(round(uptakeAc1,1))+','+str(round(uptakeAc2,1))+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle)+"\n")
  
  print("Avg.fitness(sd):\t"+str(avgfitness)+"\t"+str(sdfitness)+"\n")
  if(sdfitness>0.1):
      avgfitness=0.0
  
  return avgfitness,sdfitness
# end-def ecoliLongTermFLYCOP_oneConf
################################################################


