#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cobra #as cb
import pandas as pd
#import tabulate
import re
import sys
import getopt
import os.path
import copy
import csv
import math
import cobra.flux_analysis.variability
import subprocess
import shutil, errno
import statistics
from cobra import Reaction
import comets as c

################################################################
#loading a sbml model

def initialize_models():
    if not(os.path.exists("/home/ana/FLYCOP/MicrobialCommunities/ecoliLongTerm_TemplateOptimizeConsortiumV0/ModelsInput/iJO1366.cmd")):
        print('ERROR! Not iJO1366.cmd files with GEM of consortium strains in ModelsInput!')
    else:
        #os.chdir("ModelsInput")
        model=c.model("/home/ana/FLYCOP/MicrobialCommunities/ecoliLongTerm_TemplateOptimizeConsortiumV0/ModelsInput/iJO1366.cmd")
        print(model)
        model_id=model.id
    return model, model_id #not sure if removing model and read it in initialize_layout, just like in ecolilongtermFlYCOP...

#Stablish the layout file 
def initialize_layout(model):
    layout=c.layout(model)
    return layout
    
def ecolilongTermFLYCOP_oneConf(glu1,ac1,o21,glu2,ac2,o22, fitFunc='MaxYield_MinTime',dirPlot='',repeat=10, params_file=None):
    
    #check if class model exist
    #try:
       # var = c.model()
    #except NameError:
        #clas model doesn't exists
    model, model_id= initialize_models()
        #print(model)
    print("Fitness function:"+fitFunc)
    
  # Single GEMs parameter modifications
  # ===================================
    if not(os.path.exists('ecoli_1_tmp.cmd')):
        # 1.1.- [COBRApy] Establish modifications in model 1 
        model1=model
        model1.change_bounds('GLCtex_copy1',0,-glu1)
        model1.change_bounds('O2tex',0,-o21)
        if(ac1<=0):
            model1.change_bounds('ACtex',-1000,-ac1)
            model1.change_bounds('EX_ac_e',ac1,1000)
        else:
            model1.change_bounds('ACtex',-ac1,-1000)
            model1.change_bounds('EX_ac_e',-1000,ac1)
        model1.id='ecoli_1_tmp'
        model1.write_comets_model() #save it in the working directory ?
  
        # 1.2.- [COBRApy] Establish modifications in model 2
        model2=model
        model2.change_bounds('GLCtex_copy1',0,-glu2)
        model2.change_bounds('O2tex',0,-o22)
        if(ac1<=0):
            model2.change_bounds('ACtex',-1000,-ac2)
            model2.change_bounds('EX_ac_e',ac2,1000)
        else:
            model2.change_bounds('ACtex',-ac2,-1000)
            model2.change_bounds('EX_ac_e',-1000,ac2)
        model2.id='ecoli_2_tmp'
        model2.write_comets_model() #save it in the working directory 
    models_tmp=[model1, model2]
    layout=initialize_layout(models_tmp)
        
    
    if params_file != None:
        #establish params with a global params file
        params = c.params(global_params = params_file)
    else:
        #default params
        params = c.params()
    #establish comets class
    print(params)
    comets = c.comets(layout, params)
    comets.set_classpath('lang3', '/home/ana/comets/lib/commons-lang3-3.9/commons-lang3-3.9.jar')
    comets.set_classpath('gurobi', '/home/ana/gurobi903/linux64/lib/gurobi.jar')
    comets.set_classpath('bin', '/home/ana/comets/bin/comets_2.10.0.jar')
    #comets.JAVA_CLASSPATH= '/home/ana/Escritorio/gurobi903/linux64/lib/gurobi.jar:/home/ana/comets/lib/junit/junit-4.12.jar:/home/ana/comets/lib/junit/hamcrest-core-1.3.jar:/home/ana/comets/lib/jogl/jogamp-all-platforms/jar/jogl-all.jar:/home/ana/comets/lib/jogl/jogamp-all-platforms/jar/gluegen-rt.jar:/home/ana/comets/lib/jogl/jogamp-all-platforms/jar/gluegen.jar:/home/ana/comets/lib/jogl/jogamp-all-platforms/jar/gluegen-rt-natives-linux-amd64.jar:/home/ana/comets/lib/jogl/jogamp-all-platforms/jar/jogl-all-natives-linux-amd64.jar:/home/ana/comets/lib/JMatIO/lib/jamtio.jar:/home/ana/comets/lib/JMatIO/JMatIO-041212/lib/jmatio.jar:/home/ana/comets/lib/colt/lib/concurrent.jar:/home/ana/comets/lib/colt/lib/colt.jar:/home/ana/comets/lib/commons-lang3-3.9/commons-lang3-3.9.jar:/home/ana/comets/lib/commons-math3-3.6.1/commons-math3-3.6.1.jar:/home/ana/comets/bin/comets_2.10.0.jar'
    comets.run()
    print(comets.run_output)
      
    # To repeat X times, due to random behaviour in COMETS:
    for i in range(repeat):
        with open("output.txt", "w") as f:  #-----------don't know if needed
            comets.run()
        #graphic part
        subprocess.call(['../../Scripts/plot_biomassX2_vs_2mediaItem.sh','template', 'glc_D','ac','Ecoli1','Ecoli2'])
        
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
            biomassYieldNew=float((finalBiomass-iniBiomass)/(totGlc**0.1801559)) # molecular weigth of met1 per mmol
            # For normalizing yield
            MaximumYield=0.6
            # 7.3.- Compute second element fitnes: minimize time        
            fitTime=1-(float(endCycle)/float(240))
            # 7.4.- Compute joint fitness, as a 50% each element.
            if(fitFunc=='Yield'):
                fitness=(biomassYieldNew/MaximumYield) # Normalizing Yield
            if(float(finalBiomass-iniBiomass) > 1.03):  # Given that with both strains with WT, total biomass=1.028
                fitness=0
            
            numRxnMet=37
            pattern ="fluxes\{2\}\{1\}\{1\}\{"+i+"\}"
            with open("flux_log_template.txt", "w") as f:
                lines = f.readlines()
                for line in lines:
                    if re.search(pattern2, line):
                        line.replace(pattern2, '').replace("=", '').replace("[", '')
                        numbers=list(line.split(" "))
                        uptakeMet=numbers[numRxnMet-1]
                totalUptake.append(uptakeMet)
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

                #shutil.copy(file,file2)
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

ecolilongTermFLYCOP_oneConf(-10,-5,-18,-10,-5,-18, fitFunc='MaxYield_MinTime',dirPlot='',repeat=10, params_file=None)


# In[ ]:





# In[ ]:




