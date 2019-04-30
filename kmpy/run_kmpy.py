"""
Build and solve the ligpy kinetic model.  This model simulates the
pyrolysis of lignin in a perfectly mixed batch reactor.  For complete model
details see the journal article linked in the README at
https://github.com/houghb/ligpy.
Please read the documentation for instructions on using this module.
"""

import os
import time
import sys

import numpy as np
#import cPickle as pickle

cwd = os.getcwd()
sys.path.insert(0, cwd)
#from equivalent_compositions import write_compositionlist
#import ligpy_utils as utils
#import ddasac_utils as ddasac
from . import ode_builder

# Time the duration of running this script
#script_start_time = time.time()
#script_start_time_human = time.asctime()

# These are the files and paths that will be referenced in this program:
file_reactionlist, file_rateconstantlist, file_compositionlist\
    = ode_builder.set_paths()
working_directory = 'results_dir'
if not os.path.exists(working_directory):
    os.makedirs(working_directory)


#builiding the reactant, product and unique species list
reactants_list, products_list, unique_species = \
ode_builder.build_species_list(file_reactionlist)

#Making a complete reaction list with reactants and products
#with their stoichiometric ratios
reac_prod_list = [react + prod for react, prod in \
zip(reactants_list, products_list)]

#generating a dictionary of unique species from the species_list
#where the species are keys and the indexes are values
species_indices = {unique_species[i]:i for i in range(0, len(unique_species))}

#reversing the species_indices matrix, indexes are the keys now
indices_to_species = dict(zip(species_indices.values(), species_indices.keys()))

#building forward rate constants
temperature = 298

forward_rate_constants = ode_builder.build_kmatrix_forward(file_rateconstantlist, temperature)


print(forward_rate_constants)


