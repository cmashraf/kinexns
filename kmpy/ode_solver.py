from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import numpy as np


def initial_condition(species_list, indices_species, species_smile, concentration):
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


def stiff_ode_solver(species_list, dydtlist, y_initial, forward_rate,
                     rev_rate, t_final):
    """
    Sets up the initial condition for solving the odes
    Parameters
    ----------
    species_list     : list
                     A list of unique species in the mechanism
    dydtlist        : list
                    the list of ODEs
    y_initial       : list
                    A list of initial concentrations
    forward_rate   : list
                    A list of forward reaction rates
                    for all the reactions in the mechanism
    rev_rate        : list
                    A list of reverse reaction rates
                    for all the reactions in the mechanism
    t_final         : float
                    final time in seconds
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
    dydt = np.zeros((len(species_list)), dtype=float)
    # Define the rhs

    def rhs(t, concentration):
        kf = forward_rate
        kr = rev_rate
        y = concentration
        for i in range(len(dydtlist)):
            dydt[i] = eval(dydtlist[i].split('=')[1])
#            print(dydt)
        del t
        del y
        del kf, kr
        return dydt

    t0 = 0
    # Define an Assimulo problem
    exp_mod = Explicit_Problem(rhs, y_initial, t0)

    # Define an explicit solver
    exp_sim = CVode(exp_mod)  # Create a CVode solver

    # Sets the parameters
    exp_sim.iter = 'Newton'  # Default 'FixedPoint'
    exp_sim.discr = 'BDF'  # Default 'Adams'
    exp_sim.atol = [1e-5]  # Default 1e-6
    exp_sim.rtol = 1e-5  # Default 1e-6
    exp_sim.maxh = 0.1
#     exp_sim.display_progress = True
#     exp_sim.linear_solver = "DENSE"
    # Simulate
    t, y = exp_sim.simulate(t_final, 200)  # Simulate 5 seconds
    # return exp_mod, exp_sim
    return t, y
