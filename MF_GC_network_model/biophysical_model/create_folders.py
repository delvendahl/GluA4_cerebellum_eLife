# create folders to store temporary simulation data and results

import os
import sys


def create_folder(name):
    created = False
    if not os.path.exists(name):
        os.mkdir(name)
        created = True
    
    return created


os.chdir(sys.path[0])

# create temp folder for simulation data
foldername = 'tempdata'
create_folder(foldername)

# create folders for simulation results
foldername = 'results'
create_folder(foldername)

correlations = [0, 5, 10, 15, 20, 25, 30]

for c in correlations:
    # for original model
    foldername = 'results/orig_data_r{}'.format(c)
    create_folder(foldername)
    # for KO model
    foldername = 'results/ko_data_r{}'.format(c)
    create_folder(foldername)
