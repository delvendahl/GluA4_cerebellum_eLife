## GC integrate-and-fire model

Simulations of GC properties and MF->GC transmission from Kita et al. 2021 ([doi:10.7554/eLife.65152](https://elifesciences.org/articles/65152)).  
Simulations are run in IgorPro (Wavemetrics) using NeuroMatic (www.thinkrandom.com).  
Tested on Mac OS 10.14.6 running IgorPro 6.37 and Neuromatic v3.0c.  
  
Folder description:  
* __fixed MF input__  
Contains .ipf files to run simulations for fixed frequency MF input. These simulations run a single MF-GC synaptic input based on experimental conductances. Individual .ipf files run simulations for 100–300Hz MF input. Reproduces data from Figure 4–figure supplement 2.
* __random MF input__  
Contains .ipf files to run simulations for random (poisson) frequency MF input. These simulations run four MF-GC synaptic inputs, synaptic conductances are based on experimental results. Individual .ipf files run simulations for 10–320Hz MF input ("low_freq") or 50–1000Hz MF input ("high_freq"). Reproduces data from Figure 5 and Figure 5–figure supplement 2.  
* __tonic inhibition__  
Contains the .ipf file to run GC model with step current injections. Runs GC model with ("ctrl") or without tonic inhibition ("bmi"). Reproduces Figure 2–figure supplement 2.
* __analysis routines__  
Contains .ipf files for various analyses of simulation results.  

