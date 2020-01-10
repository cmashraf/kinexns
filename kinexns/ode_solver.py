from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import numpy as np
import warnings
import sys
import io

from assimulo.solvers.sundials import CVodeError


def initial_condition(species_list, indices_species,
                      species_smile, concentration):
    """
    Sets up the initial condition for solving the odes
    Parameters
    ----------
    species_list     : list
                     A list of unique species in the mechanism
    indices_species  : dict
                     the dictionary species_indices
    species_smile   : list
                    A list of smiles of the species with
                    initial concentrations
    concentration   : list
                    corresponding list of concentrations
                    of species_smiles
    Returns
    ----------
    matrix           : list
                    A list of initial concentrations of
                    all the species

    """
    y0 = np.zeros((len(species_list)), dtype=float)
    for i, smile in enumerate(species_smile):
        y0[indices_species[smile]] = concentration[i]
    return y0


def stiff_ode_solver(matrix, y_initial, forward_rate,
                     rev_rate, third_body=None,
                     iteration='Newton', discr='BDF',
                     atol=1e-10, rtol=1e-6,
                     sim_time=0.001, num_data_points=500):
    """
    Sets up the initial condition for solving the odes
    Parameters
    ----------
    matrix          : ndarray
                    stoichiometric matrix
    y_initial       : list
                    A list of initial concentrations
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    rev_rate        : list
                    A list of reverse reaction rates
                    for all the reactions in the mechanism
    sim_time        : float
                    total time to simulate in seconds
    third_body      : ndarray
                    third body matrix, default = None
    iteration       : str
                    determines the iteration method that is be
                    used by the solver, default='Newton'
    discr           : determines the discretization method,
                    default='BDF'
    atol            : float
                    absolute tolerance(s) that is to be used
                    by the solver, default=1e-10
    rtol            : float
                    relative tolerance that is to be
                    used by the solver, default= 1e-7
    num_data_points : integer
                    number of even space data points in output
                    arrays, default = 500
    Returns
    ----------
    t1             : list
                    A list of time-points at which
                    the system of ODEs is solved
                    [t1, t2, t3,...]
    y1              : list of lists
                    A list of concentrations of all the species
                    at t1 time-points
                    [[y1(t1), y2(t1),...], [y1(t2), y2(t2),...],...]

    """

    # y0[0] = 0
    # y0[0] = 0
    # dydt = np.zeros((len(species_list)), dtype=float)
    # Define the rhs
    kf = forward_rate
    kr = rev_rate
    mat_reac = np.abs(np.asarray(np.where(matrix < 0, matrix, 0)))
    mat_prod = np.asarray(np.where(matrix > 0, matrix, 0))

    #  d = kf * np.prod(y0**np.abs(mat_reac), axis = 1) - kr * np.prod(y0**mat_prod, axis = 1)

    def rhs(t, concentration):
        #        print(t)

        y = concentration
        if third_body is not None:
            third_body_eff = np.dot(third_body, y)
            third_body_eff = np.where(third_body_eff > 0,
                                      third_body_eff, 1)
        #            print(third_body_eff)
        else:
            third_body_eff = np.ones(len(forward_rate))
        #            print(len(third_body_matrix))
        rate_concentration = (kf * np.prod(np.power(y, mat_reac), axis=1)
                              - kr * np.prod(np.power(y, mat_prod), axis=1))

        dydt = np.dot(third_body_eff,
                      np.multiply(matrix, rate_concentration.reshape
                                  (matrix.shape[0], 1)))
        # dydt = [np.sum(third_body_eff * (matrix[:, i] * rate_concentration)) for i in
        #         range(len(species_list))]
        del t
        del y
        return dydt

    t0 = 0
    # Define an Assimulo problem
    exp_mod = Explicit_Problem(rhs, y_initial, t0)

    # Define an explicit solver
    exp_sim = CVode(exp_mod)  # Create a CVode solver

    # Sets the parameters
    exp_sim.iter = iteration  # Default 'FixedPoint'
    exp_sim.discr = discr  # Default 'Adams'
    exp_sim.atol = [atol]  # Default 1e-6
    exp_sim.rtol = rtol  # Default 1e-6
    exp_sim.maxh = 0.1
    exp_sim.minh = 1e-18
    exp_sim.num_threads = 1

    while True:
        try:
            t1, y1 = exp_sim.simulate(sim_time, num_data_points)
            break
        except CVodeError:
            # reduce absolute error by two orders of magnitude
            # and try to solve again.
            print("next process started")
            atol = atol * 1e-2
            exp_sim.atol = atol
            # if atol < 1e-15:
            #     t1 = 0
            #     y1 = 0
            #     break
            # print(exp_sim.atol)

    return t1, y1
