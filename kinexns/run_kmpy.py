"""
Build and solve the ligpy kinetic model.  This model simulates the
pyrolysis of lignin in a perfectly mixed batch reactor.  For complete model
details see the journal article linked in the README at
https://github.com/houghb/ligpy.
Please read the documentation for instructions on using this module.
"""
import io
# import time
import sys
import timeit
# from .ode_builder import *
# from .ode_solver import *
from .sensitivity_analysis import *
from .parallel_solver import *
import matplotlib.pyplot as plt
from rdkit.Chem import Descriptors
from rdkit import Chem

# import numpy as np
# from .constants import GAS_CONST, PR_ATM
# from .constants import KCAL_JL, HT_JL
# import cPickle as pickle

cwd = os.getcwd()
sys.path.insert(0, cwd)
# from equivalent_compositions import write_compositionlist
# import ligpy_utils as utils
# import ddasac_utils as ddasac

#  Time the duration of running this script
# script_start_time = time.time()
# script_start_time_human = time.asctime()

#  These are the files and paths that will be referenced in this program:
myPath = os.path.dirname(os.path.abspath(__file__))
file_reactionlist, file_rateconstantlist, file_free_energy\
    = set_paths(myPath)
working_directory = 'results_dir'
if not os.path.exists(working_directory):
    os.makedirs(working_directory)


# builiding the reactant, product and unique species list
reactants_list, products_list, unique_species =\
    build_species_list(file_reactionlist)

# Making a complete reaction list with reactants and products
# with their stoichiometric ratios
reac_prod_list = [react + prod for react, prod in
                  zip(reactants_list, products_list)]

# generating a dictionary of unique species from the species_list
# where the species are keys and the indexes are values
species_indices = {unique_species[i]: i for i in
                   range(0, len(unique_species))}

# print(species_indices)
# reversing the species_indices dictionaries, indexes are the keys now
indices_to_species = dict(zip(species_indices.values(),
                          species_indices.keys()))


# building the odes for the entire mechanism
dydt_list = build_dydt_list(reac_prod_list, unique_species,
                            species_indices, rev_rate='yes')

# building forward rate constants
temp = 773

forward_rate_constants =\
    build_kmatrix_forward(file_rateconstantlist, temp)
# print(forward_rate_constants)
# building the free energy dictionary
free_energy_dict = build_free_energy_dict(file_free_energy, temp)
# Building reverse rate constants
# gibbs_energy, change_mol = \
#   ode_builder.build_free_energy_change(reac_prod_list, free_energy_dict)
# equilibrium_constants = [np.exp(-n * 1000/(GAS_CONST * temp))
#                         for n in gibbs_energy]

# reverse_rate_constants = [(a / b) * 1000 * (GAS_CONST * temp / PR_ATM) ** c
#                          if c < 3 else 0 for (a, b, c) in
#                          zip(forward_rate_constants, equilibrium_constants,
#                          change_mol)]
reverse_rate_constants = \
    build_kmatrix_reverse(reac_prod_list, free_energy_dict,
                          forward_rate_constants, temp)

text_trap = io.StringIO()
sys.stdout = text_trap

sys.stdout = sys.__stdout__
# getting the initial conditions setup
initial_conc = initial_condition(unique_species, species_indices,
                                 ['O[C@H]1[C@H](O)CO[C@@H](O)[C@@H]1O'], [1.0])

# Solving the system of ODEs
mod, sim = stiff_ode_solver(unique_species, dydt_list, initial_conc,
                            forward_rate_constants, reverse_rate_constants)

wt_xylose = Descriptors.MolWt(Chem.MolFromSmiles('O[C@H]1[C@H](O)CO[C@@H](O)[C@@H]1O')) * initial_conc[species_indices['O[C@H]1[C@H](O)CO[C@@H](O)[C@@H]1O']]


plt.plot(mod, sim[:,species_indices['O[C@H]1[C@H](O)CO[C@@H](O)[C@@H]1O']] * 100/initial_conc[species_indices['O[C@H]1[C@H](O)CO[C@@H](O)[C@@H]1O']], color="b" , label = 'xylose')
plt.plot(mod, sim[:,species_indices['O=CCO']] * 100 * Descriptors.MolWt(Chem.MolFromSmiles('O=CCO')) / wt_xylose, color="r" , label = 'Glycolaldehyde')
plt.plot(mod, sim[:,species_indices['O=Cc1ccco1']] * 100 * Descriptors.MolWt(Chem.MolFromSmiles('O=Cc1ccco1')) / wt_xylose, color="g", label = 'Furfural')
# #P.xscale(10e-10)

plt.ylim(0, 100)
plt.legend()
plt.xlabel('Time(s)')
plt.ylabel('Wt (%)')
plt.show()
# plt.yscale('log')

# Generating parameters for Sensitivity analysis
sa_path = '/Users/chowdhury/Documents/kmpy_results/SA_data/'
gen_params(2, sa_path, 'params.txt', 'param_set.txt')

setpath = '/Users/chowdhury/Documents/kmpy_results/SA_data/test/'
file_path_read = setpath + 'param_set.txt'
file_path_write = setpath + 'model_solutions.txt'
file_input = setpath + 'params.txt'

# benchmarks = list([])
#
# benchmarks.append(timeit.Timer('serial(file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list)',
#             'from __main__ import serial, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list').timeit(number=1))
#
# benchmarks.append(timeit.Timer('multiprocess(2, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list)',
#             'from __main__ import multiprocess, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list').timeit(number=1))
#
# benchmarks.append(timeit.Timer('multiprocess(3, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list)',
#             'from __main__ import multiprocess, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list').timeit(number=1))
#
# benchmarks.append(timeit.Timer('multiprocess(4, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list)',
#             'from __main__ import multiprocess, file_path_read, forward_rate_constants, file_rateconstantlist, file_free_energy, reac_prod_list, unique_species, initial_conc, dydt_list').timeit(number=1))
#
# print(benchmarks)
