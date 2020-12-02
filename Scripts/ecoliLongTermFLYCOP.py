#!/usr/bin/python3

import cobra as cb
import pandas as pd
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
import cometspy as c
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
################################################################


def model_modifications(glu1, ac1, o21, glu2, ac2, o22):
    # Single GEMs parameter modifications
  # ===================================

    # 1.1.- [COBRApy] Establish modifications in model 1
    model1=c.model('/home/ana/FLYCOP/Scripts/iJO1366.cmd')
    model1.change_bounds('GLCtex_copy2', 0, 0)
    model1.change_bounds('PFK_3', 0, 0)
    model1.change_bounds('EX_o2_e',-1000,1000)
    model1.change_bounds('EX_glc__D_e',-1000,1000)
    model1.change_bounds('GLCtex_copy1', 0, -glu1)
    model1.change_bounds('O2tex', 0,-o21)
    if ac1<=0:
        model1.change_bounds('ACtex',-1000,-ac1)
        model1.change_bounds('EX_ac_e', ac1, 1000)
    else:
        model1.change_bounds('ACtex', -ac1,1000)
        model1.change_bounds('EX_ac_e',-1000,ac1)
    model1.id='ecoli_1_tmp'
    model1.write_comets_model()
    del(model1)

    # 1.2.- [COBRApy] Establish modifications in model 2
    model2=c.model("/home/ana/FLYCOP/Scripts/iJO1366.cmd")
    model2.change_bounds('GLCtex_copy2', 0, 0)
    model2.change_bounds('PFK_3', 0, 0)
    model2.change_bounds('EX_o2_e',-1000,1000)
    model2.change_bounds('EX_glc__D_e',-1000,1000)
    model2.change_bounds('GLCtex_copy1',0,-glu2)
    model2.change_bounds('O2tex',0,-o22)
    if(ac2<=0):
        model2.change_bounds('ACtex',-1000,-ac2)
        model2.change_bounds('EX_ac_e',ac2,1000)

    else:
        model2.change_bounds('ACtex',-ac2,1000)
        model2.change_bounds('EX_ac_e',-1000,ac2)
    model2.id='ecoli_2_tmp'
    model2.write_comets_model()
    del(model2)
    return (print('Modifications of the model done!'))

def initialize_layout(file):
    #The layout is a file with the stablished format
    layout=c.layout(file)
    return layout

def initialize_params(package, globall):
    """Function to initialize the comets' params class
    it can be initialize by submitting two files, one for the package parameters
    and one for the global ones.
    If you don't submit a file, the params class will be initialize with the values stated below
    which have been tested in this simulation"""

    if package or globall is not None:
        params = c.params(global_params = globall, package_params= package)
    else:
        params = c.params()
        params.all_params['maxCycles']=240
        params.all_params['timeStep']=0.1
        params.all_params['spaceWidth']=0.05
        params.all_params['allowCellOverlap']= True
        params.all_params['deathRate']= 0.01
        params.all_params['numRunThreads']= 4
        params.all_params['maxSpaceBiomass']= 10
        params.all_params['defaultVmax']=20
        params.all_params['showCycleTime']=True
        params.all_params['writeTotalBiomassLog']=True
        params.all_params['writeMediaLog']=True
        params.all_params['writeFluxLog']=True
        params.all_params['useLogNameTimeStamp']=False
        params.all_params['FluxLogRate']=1
        params.all_params['MediaLogRate']=1

    return params

def make_graph(comets, met, met22):

    '''This function creates a figure and saves it to pdf format.
    It also creates the file biomass_vs_met.txt which contais the quantity
    of each strain and metabolite and has the following columns:
    time(h), strain1 ... strainX, met1 ... metX.'''

    df =comets.media
    df=df.drop(columns=['x', 'y'])
    met2 =df[df['metabolite'] == str(met22)]
    met1 =df[df['metabolite'] == str(met)]
    met1 =met1.reset_index()
    met1=met1.drop(columns=['index', 'metabolite'])

    i = 45
    while i < 242:
        met1.loc[i] = [i+2, 0]
        i +=1

    met1=met1.drop(columns='cycle')



    met2=met2.reset_index()
    met2=met2.drop(columns=['metabolite', 'index'])

    j=46

    while j < 240:
        met2.loc[j] = [j+2, 0]
        j=j+1
    met2.loc[-1] = [1, 0]
    met2 =met2.sort_values(by=['cycle'])
    met2=met2.drop(columns='cycle')
    met2=met2.reset_index()
    met2=met2.drop(columns='index')
    df2=comets.total_biomass

    df2[str(met)]=met1
    df2[str(met22)]=met2

    df2.columns=['time(h)','Ecoli1', 'Ecoli2', str(met), str(met22)]
    np.savetxt(r'biomass_vs_glc_D_e_ac_e_template.txt', df2.values, fmt='%s',delimiter='\t')

    #reading the biomass file as a dataframe and adding it the metabolite concentrations
    plt.ioff()
    fig, ax = plt.subplots()
    ax.axis(ymax=1.6)
    ax.set_xlabel('time (h)')
    ax.set_ylabel('biomass (g/L)')
    l1,=ax.plot(df2['time(h)']*0.1, df2['Ecoli1'], label='Ecoli1', color='k')
    l2,=ax.plot(df2['time(h)']*0.1, df2['Ecoli2'], linestyle='--', label='Ecoli2', color='k')
    ax2 = ax.twinx()
    ax2.set_ylabel('metabolite conc (mM)')
    ax2.axis(ymin=0, ymax=10)
    l3,=ax2.plot(df2['time(h)']*0.1, df2[str(met)], label=str(met), color='b')
    l4,=ax2.plot(df2['time(h)']*0.1, df2[str(met22)], label=str(met22), color='r')
    fig.tight_layout()
    plt.legend([l1, l2, l3, l4],["Ecoli1", "Ecoli2", str(met), str(met22)])

    #Saving the graph as a pdf
    plt.savefig('biomass_vs_glc_D_e_ac_e_template_plot1.pdf')
    return df2
    #df2.drop(df2.index, inplace=True)

def end_simulation_cycle(ac1, ac2, df2, met):
    #function that stablishes the endCyle after a certain stop condition

    iniBiomass=float(df2.at[0,'Ecoli1'])+float(df2.at[0, 'Ecoli2'])
    totGlc=float(df2.at[0,'glc__D_e'])
    endGlcCycle=0
    for row in df2.iterrows():
        endCycle=int(row[1][0])
        glcConc=float(row[1][3])
        acConc=float(row[1][4])
        if endGlcCycle==0 and glcConc==0.0:
            endGlcCycle=endCycle
        if glcConc==0.0 and acConc==0.0:
            break
        if glcConc==0.0 and ac1>=0 and ac2>=0:
            break


    return iniBiomass, endCycle, totGlc

def biomass_yield(df2,endCycle):
    # Function that compute final biomass as the maximum biomass of each strain

    finalBiomass1=0
    finalBiomass2=0
    count=0
    for row in df2.iterrows():
        if float(row[1][1])>finalBiomass1:
            finalBiomass1=float(row[1][1])
        if float(row[1][2])>finalBiomass2:
            finalBiomass2=float(row[1][2])
        if count > endCycle:
            break
        count+=1
    finalBiomass=finalBiomass1+finalBiomass2
    return finalBiomass

def compute_fitness(biomassYieldNew, MaximumYield, finalBiomass, iniBiomass):
    #Function that computes the fitness according to the fitness function 'yield'
    fitness=(biomassYieldNew/MaximumYield) # Normalizing Yield
    if(float(finalBiomass-iniBiomass) > 1.03):  # Given that with both strains with WT, total biomass=1.028
        fitness=0
    return fitness

def uptakeMetabolite(comets,layout):
    #Function that calculates the uptake of the metabolite 'ac_e'
    #Changing the metabolite needs changing numRxnMet; the position of the reaction 'Ex_metName' in the model
    flux = comets.fluxes
    numRxnMet=88
    uptakeMet=[]
    i=0
    while i < len(layout.models)+1:
        uptakeMet.append(float(flux.iat[i+2,numRxnMet+3]))
        i+=1
    return uptakeMet

def ecoliLongTermFLYCOP_oneConf(glu1,ac1,o21,glu2,ac2,o22, fitFunc='MaxYield_MinTime', dirPlot='',repeat=10, params_package=None, params_global=None, layout_file='/home/ana/FLYCOP/Scripts/Nuevo_layout.txt',met='glc__D_e', met2='ac_e'):
    """Primarily function of the simulation"""

    print("Fitness function:"+fitFunc)

    if not(os.path.exists('/home/ana/FLYCOP/Scripts/ecoli_1_tmp.cmd')):
        model_modifications(glu1, ac1, o21, glu2, ac2, o22)

    layout=initialize_layout(layout_file)

    params =initialize_params(params_package, params_global)

    #establish comets class
    comets = c.comets(layout, params)

    """You may need to adjust the comets.set_classpath object if your programs aren't in the default comets location"""



    if not(os.path.exists('IndividualRunsResults')):
        os.makedirs('IndividualRunsResults')
    totfitness=0
    sumTotBiomass=0
    sumTotYield=0
    fitnessList=[]

    # To repeat X times, due to random behaviour in COMETS:
    for i in range(repeat):
	
        comets.run(delete_files=True)
	
        #obtaining the results and writing them to csv
        comets.total_biomass.to_csv('Total_biomass_log_template.txt',sep='\t',index=False)
        comets.fluxes.to_csv('Flux_log_template.txt',sep='\t', index=False)
        comets.media.to_csv('Media_log_template.txt',sep='\t', index=False)

        #graphic part
        df2 =make_graph(comets, met, met2)

        # 7.- Compute fitness (measure to optimize):
        print('computing fitness...')

        # 7.1.- Determine endCycle: when glucose and acetate are exhausted
        iniBiomass, endCycle, totGlc = end_simulation_cycle(ac1, ac2, df2, met)

        # 7.2.- Compute first element fitness: maximize biomass yield
        finalBiomass = biomass_yield(df2,endCycle)

        biomassYieldNew=float((finalBiomass-iniBiomass)/(totGlc*0.1801559)) # molecular weigth of met1 per mmol

        # For normalizing yield
        MaximumYield=0.6

        # 7.3.- Compute second element fitnes: minimize time
        fitTime=1-(float(endCycle)/float(240))

        # 7.4.- Compute joint fitness, as a 50% each element.
        fitness= compute_fitness(biomassYieldNew, MaximumYield, finalBiomass, iniBiomass)

        # 7.5.- Calculate the uptake of ac_e
        uptakeMet = uptakeMetabolite(comets,layout)

        print(" Total biomass: "+str(round(finalBiomass,6))+" in cycle "+str(endCycle)+". Biomass yield="+str(round(biomassYieldNew,6)))
        totfitness=totfitness+fitness
        fitnessList.append(fitness)
        sumTotBiomass=sumTotBiomass+finalBiomass
        sumTotYield=sumTotYield+biomassYieldNew

        # Copy individual solution
        file='IndividualRunsResults/'+'biomass_vs_glc_D_ac_run'+str(i)+'_'+str(fitness)+'_'+str(endCycle)+'.pdf'
        shutil.move('biomass_vs_glc_D_e_ac_e_template_plot1.pdf',file)
        if(dirPlot != ''):
            file2=dirPlot+'biomass_vs_glc_D_ac_'+str(glu1)+'_'+str(ac1)+'_'+str(o21)+'_'+str(glu2)+'_'+str(ac2)+'_'+str(o22)+'_'+str(round(uptakeMet[0],1))+'_'+str(round(uptakeMet[1],1))+'_run'+str(i)+'_'+str(fitness)+'_'+str(endCycle)+'.pdf'
            shutil.copy(file,file2)



        file='IndividualRunsResults/'+'total_biomass_log_run'+str(i)+'.txt'
        shutil.move('Total_biomass_log_template.txt',file)
        file='IndividualRunsResults/'+'media_log_run'+str(i)+'.txt'
        shutil.move('Media_log_template.txt',file)
        file='IndividualRunsResults/'+'flux_log_run'+str(i)+'.txt'
        shutil.move('Flux_log_template.txt',file)

    avgfitness=totfitness/repeat
    sdfitness=statistics.stdev(fitnessList)
    avgBiomass=sumTotBiomass/repeat
    avgYield=sumTotYield/repeat
    print("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle")
    print(fitFunc+"\t"+str(glu1)+','+str(ac1)+','+str(o21)+','+str(glu2)+','+str(ac2)+','+str(o22)+','+str(round(uptakeMet[0],1))+','+str(round(uptakeMet[1],1))+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle))
    with open(dirPlot+"configurationsResults"+fitFunc+".txt", "a") as myfile:
        myfile.write("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle\n")
        myfile.write(fitFunc+"\t"+str(glu1)+','+str(ac1)+','+str(o21)+','+str(glu2)+','+str(ac2)+','+str(o22)+','+str(round(uptakeMet[0],1))+','+str(round(uptakeMet[1],1))+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle)+"\n")
    print("Avg.fitness(sd):\t"+str(avgfitness)+"\t"+str(sdfitness)+"\n")
    if(sdfitness>0.1):
        avgfitness=0.0
    try:
        os.remove('/home/ana/FLYCOP/Scripts/ecoli_1_tmp.cmd')
    except FileNotFoundError:
        pass
    try:
        os.remove('/home/ana/FLYCOP/Scripts/ecoli_2_tmp.cmd')
    except FileNotFoundError:
        pass
    return avgfitness,sdfitness
