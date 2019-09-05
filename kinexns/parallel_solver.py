import multiprocessing as mp
from .ode_builder import *
from .ode_solver import *
from .parse_chemkin import *
import numpy as np
import os
import pandas as pd


def func_solv(data, forward_rate, file_rateconstant, file_energy,
              matrix, species_list, initial_y, t_final, factor,
              third_body, pos=None, chemkin_data=None,
              smiles=None):
    """
    Solves the system of ODEs for different rate constants
    generated from the data file
    Parameters
    ----------
    data             : str
                    Each line of the data files where various randomly chosen
                    parameter values are listed to perform sensitivity analysis
    matrix          : ndarray
                    stoichiometric matrix
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    file_rateconstant : str
                      path to the file `complete_rateconstantlist.dat`
    file_energy     : str
                    path to the file 'free_energy_library.dat'
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    smiles             : dict
                    the smiles dictionary generated from
                    species_smiles.dat file
    t_final        : float
                    final time in seconds
    factor              : float
                        conversion factor from given unit of energy to kJ
    third_body          : ndarray
                        matrix with third body efficiencies
    pos                 : int
                        position argument for multiprocessing
    chemkin_data        :ndarray
                        the data from parsed chemkin reaction file
    Returns
    ----------
    sim[-1]         : list
                    A list of final concentrations of all the
                    species at t_final
    """
    kf_random = np.zeros((len(forward_rate)), dtype=float)
    values = list(map(float, data.split()))
    kf_random[:] = np.array(values[:len(forward_rate)])
    # temp = values[-1] + 273
    temp = 950.0
    kf_actual = np.array(build_forward_rates(file_rateconstant, temp))
    kf_pur = np.array([actual * 10 ** rand for actual,
                       rand in zip(kf_actual, kf_random)])
    kf_purturbed = list(kf_pur)
    if chemkin_data:
        chemkin = True
        free_energy_dict = generate_thermo_dict(file_energy, smiles, temp)
        kf_purturbed = \
            update_rate_constants_for_pressure(chemkin_data,
                                               kf_purturbed, temp)
    else:
        chemkin = False
        free_energy_dict = build_free_energy_dict(file_energy, temp)
    kr_purturbed = build_reverse_rates(free_energy_dict, species_list, matrix,
                                       factor, kf_purturbed, temp, chemkin)
    mod, sim = stiff_ode_solver(matrix, initial_y, kf_purturbed,
                                kr_purturbed, third_body=third_body,
                                sim_time=t_final)

    return pos, sim[-1]


def serial_ss(file_read, forward_rate, file_rateconstant,
              file_energy, matrix, species_list, factor,
              initial_y, t_final, third_body=None,
              chemkin_data=None, smiles=None, chemkin=True):
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
    matrix          : ndarray
                    stoichiometric matrix
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    t_final         : float
                    final time in seconds
    third_body          : ndarray
                        matrix with third body efficiencies
    chemkin_data        :ndarray
                        the data from parsed chemkin reaction file
    smiles             : dict
                    the smiles dictionary generated from
                    species_smiles.dat file
    factor              : float
                        conversion factor from given unit of energy to kJ
    chemkin             : bool
                        indicates if chemkin files are read as input files
                        default = True
    Returns
    ----------
                     : list
                    A list of final concentrations of all the
                    species at t_final for all the given combinations
                    of parameters listed in 'param_set.txt' file
    """
    read_file = open(file_read, "r")
    results = []
    for pos, data in enumerate(read_file):
        result = func_solv(data, forward_rate, file_rateconstant, file_energy,
                           matrix, species_list, initial_y, t_final, factor,
                           third_body, pos, chemkin_data, smiles)
        results.append(result)
    return results


def multiprocess(processes, file_read, forward_rate, file_rateconstant,
                 file_energy, matrix, species_list, factor,
                 initial_y, t_final, third_body=None,
                 chemkin_data=None, smiles=None, chemkin=True):
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
    matrix          : ndarray
                    stoichiometric matrix
    file_energy     : str
                    path to the file 'free_energy_library.dat'
    species_list     : list
                     A list of unique species in the mechanism
    initial_y       : list
                    A list of initial concentrations
    t_final         : float
                    final time in seconds
    third_body          : ndarray
                        matrix with third body efficiencies
    chemkin_data        :ndarray
                        the data from parsed chemkin reaction file
    smiles             : dict
                    the smiles dictionary generated from
                    species_smiles.dat file
    factor              : float
                        conversion factor from given unit of energy to kJ
    chemkin             : bool
                        indicates if chemkin files are read as input files
                        default = True
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
        matrix, species_list, initial_y, t_final, factor,
        third_body, pos, chemkin_data, smiles))
               for (pos, data) in enumerate(read_file)]
    results = [p.get() for p in results]
    results.sort()  # to sort the results by input window width
    results = [r[1] for r in results]
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
    dataframe.to_csv(file_name, sep='\t', header=None, index=None)
