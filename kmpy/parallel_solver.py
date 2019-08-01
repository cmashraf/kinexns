import multiprocessing as mp
from .ode_builder import *
from .ode_solver import *
import numpy as np
import os


def func_solv(data, forward_rate, file_rateconstant, file_energy,
              complete_list, species_list, initial_y, dydtlist, t_final):
    """
    Solves the system of ODEs for different rate constants
    generated from the data file
    Parameters
    ----------
    data             : str
                    Each line of the data files where various randomly chosen
                    parameter values are listed to perform sensitivity analysis
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    file_rateconstant : str
                      path to the file `complete_rateconstantlist.dat`
    file_energy     : str
                    path to the file 'free_energy_library.dat'
    complete_list   : list
                     A list of all the reactions with reactant and
                      product species and their stoichimetric coeffs
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    dydtlist        : list
                    the list of ODEs
    t_final         : float
                    final time in seconds
    Returns
    ----------
    sim[-1]         : list
                    A list of final concentrations of all the
                    species at t_final
    """
    kf_random = np.zeros((len(forward_rate)), dtype=float)
    values = list(map(float, data.split()))
    kf_random[:] = np.array(values[:len(forward_rate)])
    temp = values[-1] + 273
    kf_actual = np.array(build_kmatrix_forward(file_rateconstant, temp))
    kf_pur = np.array([actual * 10 ** rand for actual,
                      rand in zip(kf_actual, kf_random)])
    kf_purturbed = kf_pur.tolist()
    free_energy_dict = build_free_energy_dict(file_energy, temp)
    kr_purturbed = build_kmatrix_reverse(complete_list, free_energy_dict, kf_purturbed, temp)
    mod, sim = stiff_ode_solver(species_list, dydtlist, initial_y,
                                kf_purturbed, kr_purturbed, t_final)

    return sim[-1]


def serial(file_read, forward_rate, file_rateconstant, file_energy,
           complete_list, species_list,
           initial_y,  dydtlist, t_final):
    """
    Iteratively solves the system of ODEs for different rate constants
    generated from the data file in serial
    Parameters
    ----------
    file_read       : str
                    path of the 'param_set' file where all
                    the parameter combinations are listed
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    file_rateconstant : str
                      path to the file `complete_rateconstantlist.dat`
    file_energy     : str
                    path to the file 'free_energy_library.dat'
    complete_list   : list
                     A list of all the reactions with reactant and
                      product species and their stoichimetric coeffs
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    dydtlist        : list
                    the list of ODEs
    t_final         : float
                    final time in seconds
    Returns
    ----------
                     : list
                    A list of final concentrations of all the
                    species at t_final for all the given combinations
                    of parameters listed in 'param_set.txt' file
    """
    read_file = open(file_read, "r")
    return [func_solv(data, forward_rate, file_rateconstant, file_energy,
                      complete_list, species_list,
                      initial_y, dydtlist, t_final) for data in read_file]


def multiprocess(processes, file_read, forward_rate, file_rateconstant,
                 file_energy, complete_list,
                 species_list, initial_y, dydtlist, t_final):

    """
    Iteratively solves the system of ODEs for different rate constants
    generated from the data file in parellel
    Parameters
    ----------
    processes       : int
                    number of processors
    file_read       : str
                    path of the 'param_set' file where all
                    the parameter combinations are listed
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    file_rateconstant : str
                      path to the file `complete_rateconstantlist.dat`
    file_energy     : str
                    path to the file 'free_energy_library.dat'
    complete_list   : list
                     A list of all the reactions with reactant and
                      product species and their stoichimetric coeffs
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    dydtlist        : list
                    the list of ODEs
    t_final         : float
                    final time in seconds
    Returns
    ----------
    results         : list
                    A list of final concentrations of all the
                    species at t_final for all the given combinations
                    of parameters listed in 'param_set.txt' file
    """
    read_file = open(file_read, "r")
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(func_solv, args=(
               data, forward_rate, file_rateconstant, file_energy,
               complete_list, species_list, initial_y,
               dydtlist, t_final)) for data in read_file]
    results = [p.get() for p in results]
    results.sort()  # to sort the results by input window width
    return results


def write_model_sol(file_name, res):
    """
    Writes the model solutions in a file
     Parameters
    ----------
    file_name       : str
                    name of the file to save with path appended
    res             : list
                    A list of model solutions generated using
                    either serial or parallel solver
    """
    try:
        os.remove(file_name)
    except OSError:
        pass

    dataframe = pd.DataFrame.from_records(res)
    dataframe.to_csv(file, sep='\t', header=None, index=None)
