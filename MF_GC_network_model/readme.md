## Feedforward MF-GC model

Computational model of a feedforward MF-GC network. Original model taken from Cayco-Gajic et al. 2017 (doi:10.1038/s41467-017-01109-y). Implemented to run simulations of MF-GC network in Kita et al. 2021 (doi:xxx). The model simulates 640 input patterns of MF activity and the output of GCs in response to these patterns. Simualtions are performed using the parameters of the original model ('original'), and using scaled conductances based on experimental data obtained from GluA4-KO mice ('KO').  
  
Simulations were run in Python 2.7 using pyNeuroML 0.5.5 and jNeuroML 0.8.5.  
Dependencies:  
* pyNeuroML (https://github.com/NeuroML/pyNeuroML)  
* jNeuroML (https://github.com/NeuroML/jNeuroML)  
* numpy, matplotlib, scipy  
  
Model definitions, network structure and analysis routines are taken from Cayco-Gajic et al. 2017 with minor modifications. Runs 640 input patterns per condition.

Folder description:  
* __biophysical_model__  
Contains python files required to run the simulation and analyse data. Simulation output is saved into subfolders.
* __grc_lemsDefinitions__  
Contains the LEMS definition files for the integrate-and-fire GC and the AMPA and NMDA synapses.
* __input_statistics__  
Contains the MF input patterns for differenc levels of spatial correlation. Files for correlation radii of 5,15,20,25,30 are taken from Cayco-Gajic et al. 2017.
* __network_structures__  
Contains MF-GC connectivity files. Note that this simulation only uses a fixed MF-GC connectivity ratio of 4 (i.e. n_syn=4 of Cayco-Gajic et al. 2017). 'plot_network.py' can be used to generate a plot of the network structure.  
  
Getting started:  
* install dependencies  
* create folders to store data within the biophysical_model folder  

Run simulations:  
* `cd` into biophysical_model folder  
* initialize network by running `initialize_network.py`; prints total number of runs  
* run the simulation as `run_mf_gc_network.py`  
* data will be saved into subfolders named data_r0, data_r5, etc.  
  
Analyse data:  
* extract spike times from .dat files by running `save_samples_as_txt.py foldername` (replace `foldername` according to target folder); this creates .txt files in the target folder   
* analyse population characteristics by running `get_spar_cov.py foldername`(replace `foldername` according to target folder); generates file xxx.txt in the target folder  
* analyse learning performance by running `run_learning.py foldername`(replace `foldername` according to target folder); generates file xxx.txt in the target folder  
