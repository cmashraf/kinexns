import multiprocessing as mp
from .ode_builder import *
from .ode_solver import *
import numpy as np


def func_solv_new(data, forward_rate, file_rateconstant, file_energy,
                  complete_list, uniquespecies, initial_y, dydtlist):

    kf_random = np.zeros((len(forward_rate)), dtype=float)
    values = list(map(float, data.split()))
    kf_random[:] = np.array(values[:len(forward_rate)])
    temp = values[-1]
    kf_actual = np.array(build_kmatrix_forward(file_rateconstant, temp))
    kf_purturbed = np.array([actual * 10 ** rand for actual, rand in zip(kf_actual, kf_random)])
    free_energy_dict = build_free_energy_dict(file_energy, temp)
    kr_purturbed = build_kmatrix_reverse(complete_list, free_energy_dict, kf_purturbed, temp)
    mod, sim = stiff_ode_solver(uniquespecies, dydtlist, initial_y,
                                kf_purturbed, kr_purturbed)

    return sim[-1]


def serial(file_read, forward_rate, file_rateconstant, file_energy,
           complete_list, uniquespecies, species_indices, dydtlist):

    read_file = open(file_read, "r")
    return [func_solv_new(data, forward_rate, file_rateconstant, file_energy,
                          complete_list, uniquespecies,
                          species_indices, dydtlist) for data in read_file]


def multiprocess(processes, file_read, forward_rate, file_rateconstant,
                 file_energy, complete_list,
                 uniquespecies, species_indices, dydtlist):

    read_file = open(file_read, "r")
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(func_solv_new, args=(
               data, forward_rate, file_rateconstant, file_energy, complete_list,
               uniquespecies, species_indices,
               dydtlist)) for data in read_file]
    results = [p.get() for p in results]
    results.sort(key=lambda x: x[0])  # to sort the results by input window width
    return results
