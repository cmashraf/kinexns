import numpy as np
import spotpy
from .ode_solver import *


class SpotpySetup(object):

    def __init__(self, opt_parameters, sp_indices, opt_dist, cost_func,
                 forward_rate, file_rateconstant, file_energy, matrix,
                 species_list, initial_y, t_final, factor, third_body,
                 algorithm, opt_species='all', pos=None,
                 chemkin_data=None, smiles=None):
        self.parameternames = opt_parameters
        self.speciesnames = opt_species
        self.cost_func = cost_func
        self.dim = len(self.parameternames)
        self.params = []
        self.forward_rates = forward_rate
        self.opt_dist = opt_dist
        self.sp_indices = sp_indices
        self.rateconstant = file_rateconstant
        self.energy = file_energy
        self.matrix = matrix
        self.species_list = species_list
        self.initial_y = initial_y
        self.t_final = t_final
        self.algorithm = algorithm
        self.factor = factor
        self.third_body = third_body
        self.pos = pos
        self.chemkin_data = chemkin_data
        self.smiles = smiles

        for i in range(self.dim):
            self.params.append(spotpy.parameter.Uniform(self.parameternames[i], 0.1, 2, 0.2, 0.5, 0.1, 2))

        self.obs = []
        if self.speciesnames == 'all':
            self.obs = self.opt_dist
        else:
            for sp in self.speciesnames:
                self.obs.append(self.opt_dist[:, self.sp_indices[sp]])
        self.opt_obs = np.array(self.obs).flatten()

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        x = np.array(vector)
        simulations1 = np.array(stiff_ode_solver(x))

        simulations = [i / (0.1 * j) if j > 0.0 else i for i, j in zip(simulations1, observations1)]
        #        print(simulations)
        return simulations

    def evaluation(self):
        #   observations = np.array(solve_ode_actual())
        observations = [j / (0.1 * j) if j > 0.0 else 0.0 for j in self.opt_obs]
        return observations

    def objectivefunction(self, simulation, evaluation):
        if self.algorithm == 'abc' or self.algorithm == 'fscabc':
            objectivefunction = getattr(spotpy.objectivefunctions,
                                        self.cost_func)(evaluation,
                                                        simulation)
        else:
            objectivefunction = - getattr(spotpy.objectivefunctions,
                                          self.cost_func)(evaluation,
                                                          simulation)
        return objectivefunction
