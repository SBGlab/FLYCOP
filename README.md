---
Version: "FLYCOP 1.5"
author:
  - "Ana del Ramo Galian"
  - "David San Leon Granado"
  - "Beatriz Garcia-Jimenez"
---

![](FLYCOP_logo.jpg)

# FLYCOP 1.5

FLYCOP (FLexible sYnthetic Consortium OPtimization) is a framework that improves the understanding of the metabolic behaviour of microbial consortia and to automatize the modeling of those communities, by designing and optimizing enginered microbial consortia given a particular goal.

FLYCOP contributes with multiple and assorted applications, such as simulating different scenarios before in-vivo experiments; defining medium composition and detecting limiting nutrients; discovering the biological metric optimized in an evolutionary process; optimizing cross-feeding relationships; optimizing strain ratios in the consortium; etc.

**Citation**: This repository contains the code and configuration files reproducing the study cases described in (please, cite us if you use FLYCOP in your work):

**Beatriz García-Jiménez, José Luis García, Juan Nogales; FLYCOP: metabolic modeling-based analysis and engineering microbial communities, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i954–i963, [doi: 10.1093/bioinformatics/bty561](https://doi.org/10.1093/bioinformatics/bty561)**


***
### Installation

#### (a) Your-self installation: basic pre-requisites

FLYCOP software run in LINUX OS. FLYCOP can be run installing the pre-requisites individual software by yourself. 
Define the location of your personal [gurobi](http://www.gurobi.com/academia/for-universities) solver license (required by COMETS) in the container (for example, \<path_to_gurobi_license\>=/home/user):
```{sh eval=FALSE}
GRB_LICENSE_FILE=<path_to_gurobi_license>/gurobi.lic
```

FLYCOP pipeline uses some software (and all their dependencies), which must be installed before:

* [COMETS](http://www.bu.edu/segrelab/comets/) (v2.10) (faster with gurobi solver)  

Additionally, [R software](https://www.r-project.org/) is required.

#### Installation
```
#Create a conda environment
conda create -n FLYCOP python=3.8 pip
# Activate the environment
conda activate FLYCOP
# Install SWIG requirement
conda install gxx_linux-64 gcc_linux-64 swig
# INSTALL FLYCOP package
pip install FLYCOP
```
***
### Input and output description

After required software installation, you can run FLYCOP for the microbial consortia where FLYCOP was applied with the files provided in this site. And you can check the template for the FLYCOP pipeline provided to modify with an specific microbial consortium.
If you want to apply FLYCOP to design and optimize a new consortium, you could take as template the available files from one of the optimized consortia, and then to code the following required inputs to FLYCOP, described here:

####  *INPUT:*
  
  1. **FLYCOP pipeline _ad-hoc_ for a specific microbial consortium** [consortiumPrefixFLYCOP.py]: it means a python script to run a single configuration, including how to: 1) dynamically update single models and community parameters depending on the parametrized consortium configuration, 2) simulate that configured consortium in a dynamic way, and 3) evaluate the quality of the given consortium using a fitness function. This file also must include the method *initialize_models()* to update the original single strain genome scale metabolic models (GEMs) to use as base in the consortium optimization.
  2. **ConsortiumPrefix_TemplateFolder** [ConsortiumPrefix_TemplateOptimizeConsortiumV\<A>]: a directory including a layout file with the culture medium definition, the original GEMs (in matlab format) and the consortium simulation configuration.
  3. **Optimization configuration** [consortiumPrefix_confFLYCOP_scenario_v\<Y>.txt]: it defines the number of consortium configurations to evaluate (numberOfRunsLimit) and identifies the two files describing the consortium optimization:  
    + **Parameter values** [consortiumPrefix_confFLYCOP_params_v\<Z>.pcs]: it lists the range of values per parameter (in SMAC format), among those to choose by the Stochastic Local Search procedure.  
    + **Wrapper file** [consortiumPrefix_wrapperFLYCOP_v\<Y>.py]: to select the fitness function to optimize, to define the 'version identifier' for all results files, and to call the FLYCOP pipeline.  
  
"ConsortiumPrefix_TemplateFolder" must be located in *MicrobialCommunities/* directory, and the remaining input files in *Scripts/* directory.

##### Genome-scale models

GEMs used by FLYCOP cases of study can be obtained from [BiGG models database](http://bigg.ucsd.edu/) or from their respective publications (in SBML format):  

* iJO1366 [[Orth et al.,2011]](https://doi.org/10.1038/msb.2011.65)  

####  *OUTPUT:*  
FLYCOP provides different resources for robustness, sensitivity and data analysis support, being the most relevant the following ones:  

* Best configuration given the strains, media, fitness function and parameter configuration  
* Scatterplot showing explored values by each parameter  
* Correlation values and ellipse plots between different parameter and fitness values 
* Tab file with all configurations including parameter and fitness values, and some other interesting metrics (such as medium concentration of some relevant metabolites). This output would be important for further data analysis.  
* Growth curves of all explored consortium configurations  

***


