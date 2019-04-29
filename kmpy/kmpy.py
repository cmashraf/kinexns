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

#from equivalent_compositions import write_compositionlist
#import ligpy_utils as utils
#import ddasac_utils as ddasac
import ode_builder

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
reactant_list, product_list, unique_species = \
ode_builder.build_species_list(file_reactionlist)

print(unique_species)
