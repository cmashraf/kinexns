import numpy as np
import spotpy
import multiprocessing as mp
from .ode_solver import *
from .ode_builder import *
from .parse_chemkin import *


class SpotpySetup(object):

    def __init__(self, opt_params, initial_val, sp_indices,
                 opt_dist, cost_function, forward_rate, rateconstant_file,
                 energy_file, matrix, species_list, initial_y,
                 final_t, third_body, algorithm, temper, opt_species='all',
                 chemkin_data=None, smiles=None, factor=0.001):

        self.parameternames = opt_params
        self.initial_guess = initial_val
        self.speciesnames = opt_species
        self.cost_func = cost_function
        self.dim = len(self.parameternames)
        self.params = []
        self.forward_rates = forward_rate
        self.opt_dist = opt_dist
        self.sp_indices = sp_indices
        self.rateconstant = rateconstant_file
        self.energy = energy_file
        self.matrix = matrix
        self.species_list = species_list
        self.initial_y = initial_y
        self.t_final = final_t
        self.algorithm = algorithm
        self.factor = factor
        self.third_body = third_body
        self.temp = temper
        # self.pos = pos
        self.chemkin_data = chemkin_data
        self.smiles = smiles
        self.factor = factor
        self.test_species = [item for item in self.species_list if
                             item not in self.speciesnames]

        # generate parameter sets for the given parameters
        for i in range(self.dim):
            self.params.append(spotpy.parameter.Uniform
                               (self.parameternames[i], -1, 1, 0.02, 0.0))

        self.obs = []
        self.test_obs = []
        if self.speciesnames == 'all':
            self.obs = np.array(self.opt_dist)
            self.test_obs = np.array([0])
        else:
            self.opt_dist = np.array(self.opt_dist)
            for sp in self.speciesnames:
                self.obs.append(self.opt_dist[:, self.sp_indices[sp]])
            for sp in self.test_species:
                self.test_obs.append(self.opt_dist[:, self.sp_indices[sp]])
        self.opt_obs = np.array(self.obs).flatten()
        self.opt_test_obs = np.array(self.test_obs).flatten()

        self.database = open('{}.txt'.format(self.algorithm), 'w')

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # initiate the forward rate constants with the initial values
        x = np.array(vector)
        forward_rate_param = self.forward_rates.copy()
        for i, parm in enumerate(self.parameternames):
            # print(len(x))

            # print(forward_rate_param[230])
            forward_rate_param[int(parm[1:])] = \
                np.multiply(np.power(10, x[i]), forward_rate_param[int(parm[1:])])
            # print(forward_rate_param[230])
        if self.chemkin_data:
            chemkin = True
            free_energy_dict = generate_thermo_dict(self.energy,
                                                    self.smiles, self.temp)
            forward_rate_param = \
                update_rate_constants_for_pressure(self.chemkin_data,
                                                   forward_rate_param,
                                                   self.temp)
        else:
            chemkin = False
            free_energy_dict = build_free_energy_dict(self.energy, self.temp)

        rev_rate = build_reverse_rates(free_energy_dict, self.species_list,
                                       self.matrix, self.factor,
                                       forward_rate_param, self.temp, chemkin)

        time, simulations1 = stiff_ode_solver(self.matrix, self.initial_y,
                                              forward_rate_param, rev_rate,
                                              self.third_body,
                                              sim_time=self.t_final)

        results = []
        test_results = []
        if self.speciesnames == 'all':
            results = np.array(simulations1)
            test_results = np.array([0])
        else:
            simulations1 = np.array(simulations1)
            for sp in self.speciesnames:
                results.append(simulations1[:, self.sp_indices[sp]])
            for sp in self.test_species:
                test_results.append(simulations1[:, self.sp_indices[sp]])
        results = np.array(results).flatten()
        test_results = np.array(test_results).flatten()
        # print(results.shape)
        simulations = [i / (0.1 * j) if j > 0.0 else i for i, j in
                       zip(results, self.opt_obs)]
        test_simulations = [i / (0.1 * j) if j > 0.0 else i for i, j in
                            zip(test_results, self.opt_obs)]
        #        print(simulations)
        return simulations, test_simulations

    def evaluation(self):
        #   observations = np.array(solve_ode_actual())
        observations = [j / (0.1 * j) if j > 0.0 else 0.0
                        for j in self.opt_obs]
        test_observations = [j / (0.1 * j) if j > 0.0 else 0.0
                             for j in self.opt_test_obs]
        return observations

    def objectivefunction(self, simulation, evaluation):
        #         print(self.algorithm)
        if self.cost_func == 'sae':
            if self.algorithm in ['abc', 'fscabc']:
                objectivefunction = sae_func(evaluation[0], simulation[0])
                test_objectivefunction = sae_func(evaluation[1],
                                                  simulation[1])
            else:
                objectivefunction = - sae_func(evaluation[0], simulation[0])
                test_objectivefunction = - sae_func(evaluation[1],
                                                    simulation[1])

        else:
            if self.algorithm in ['abc', 'fscabc']:
                objectivefunction = getattr(spotpy.objectivefunctions,
                                            self.cost_func)(evaluation[0],
                                                            simulation[0])
                test_objectivefunction = \
                    getattr(spotpy.objectivefunctions,
                            self.cost_func)(evaluation[1], simulation[1])
            else:
                objectivefunction = - getattr(spotpy.objectivefunctions,
                                              self.cost_func)(evaluation[0],
                                                              simulation[0])
                test_objectivefunction = \
                    - getattr(spotpy.objectivefunctions,
                              self.cost_func)(evaluation[1], simulation[1])

        return objectivefunction

    def save(self, objectivefunctions, parameter, simulations, chains=None):
        parameter = list(parameter)
        line = str(objectivefunctions) + ',' + str(parameter).strip('[]') + '\n'
        self.database.write(line)


def sae_func(predictions, targets):
    return ((np.array(predictions) - np.array(targets)) ** 2).sum()


def optimization(pos, rep, opt_params, initial_val, sp_indices,
                 opt_dist, cost_function, forward_rate, rate_file,
                 energy_file, matrix, species_list, initial_y,
                 final_t, third_body, algorithm, temper, opt_species='all',
                 factor=0.001, chemkin_data=None, smiles=None):
    print(algorithm)
    parallel = "seq"
    dbformat = "custom"
    timeout = 1e6
    spot_setup = SpotpySetup(opt_params, initial_val, sp_indices, opt_dist, cost_function,
                             forward_rate, rate_file, energy_file, matrix,
                             species_list, initial_y, final_t, third_body,
                             algorithm, temper, opt_species,
                             factor=factor, chemkin_data=chemkin_data, smiles=smiles)

    sampler = getattr(spotpy.algorithms, algorithm)(spot_setup,
                                                    dbformat=dbformat)

    # print(sampler)
    sampler.sample(rep)
    spot_setup.database.close()
    # if algorithm == 'sa':
    #    sampler.sample(repetation, Tini=200)
    result = sampler.getdata()
    # print(results)
    return pos, result


def multi_optimization(processes, rep, opt_params, initial_val,
                       sp_indices, opt_dist, cost_function, forward_rate,
                       rate_file, energy_file, matrix, species_list, initial_y,
                       final_t, algorithms, temper, opt_species='all', third_body=None,
                       chemkin_data=None, smiles=None, factor=0.001):
    print(processes)
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async
               (optimization, args=(pos, rep, opt_params, initial_val, sp_indices,
                                    opt_dist, cost_function, forward_rate, rate_file,
                                    energy_file, matrix, species_list, initial_y,
                                    final_t, third_body, al, temper, opt_species,
                                    factor, chemkin_data, smiles))
               for (pos, al) in enumerate(algorithms)]
    results = [p.get() for p in results]
    results.sort()  # to sort the results by input window width
    results = [r[1] for r in results]
    return results
