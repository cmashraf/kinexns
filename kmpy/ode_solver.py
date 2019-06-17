from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import numpy as np


def initial_condition(uniquespecies, speciesindices, species_smile, concentration):

    y0 = np.zeros((len(uniquespecies)), dtype=float)
    for i, smile in enumerate(species_smile):
        y0[speciesindices[smile]] = concentration[i]
    return y0


def stiff_ode_solver(specieslist, dydtlist, y_initial, forward_rate, rev_rate):
    r"""
    Demonstration of the use of CVode by solving the
    linear test equation :math:`\dot y = - y`

    on return:

       - :dfn:`exp_mod`    problem instance

       - :dfn:`exp_sim`    solver instance

    """

    # y0[0] = 0
    dydt = np.zeros((len(specieslist)), dtype=float)
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
    t1, y1 = exp_sim.simulate(10)  # Simulate 5 seconds
    # return exp_mod, exp_sim
    return t1, y1
