import numpy as np
import matplotlib.pyplot as plt
import copy

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

