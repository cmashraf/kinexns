import numpy as np
import spotpy
import multiprocessing as mp
from .ode_solver import *
from .ode_builder import *
from .parse_chemkin import *
from .constants import KCAL_JL, HT_JL, CAL_JL
from .constants import GAS_CONST, PR_ATM


class ParamOptimize(object):

    def __init__(self, reaction_list, opt_type, sp_indices, ini_val,
                 opt_dist, cost_function, forward_rate, rateconstant_file,
                 energy_file, matrix, species_list, initial_y,
                 final_t, third_body, algorithm, temper, conver, opt_species='all',
                 chemkin_data=None, smiles=None, factor=1, pos=1):

        self.reac_list = reaction_list
        self.parameternames = create_parameter_name_list(self.reac_list, opt_type)
        self.speciesnames = opt_species
        self.cost_func = cost_function
        self.dim = len(self.parameternames)
        self.params = []
        self.opt_type = opt_type
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
        self.conv = conver
        self.test_species = [item for item in self.species_list if
                             item not in self.speciesnames]
        file_name = str(pos) + '_' + '{}.csv'.format(self.algorithm)
        self.database = open(file_name, 'w')

        # generate parameter sets for the given parameters
        for i in range(len(self.reac_list)):
            if self.opt_type == 'params':
                self.params.append(spotpy.parameter.Uniform
                                   (self.parameternames[3 * i], ini_val[0], ini_val[1], ini_val[2], ini_val[3]))
                self.params.append(spotpy.parameter.Uniform
                                   (self.parameternames[3 * i + 1], ini_val[0], ini_val[1], ini_val[2], ini_val[3]))
                self.params.append(spotpy.parameter.Uniform
                                   (self.parameternames[3 * i + 2], ini_val[0], ini_val[1], ini_val[2], ini_val[3]))
            else:
                self.params.append(spotpy.parameter.Uniform
                                   (self.parameternames[i], ini_val[0], ini_val[1], ini_val[2], ini_val[3]))
        self.opt_obs, self.opt_test_obs = get_flatten_data(self.speciesnames, self.test_species, self.opt_dist,
                                                           self.temp, self.sp_indices)

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # initiate the forward rate constants with the initial values
        x = np.array(vector)

        forward_rate_val = self.forward_rates[0].copy()
        forward_rate_params = self.forward_rates[1].copy()
        # forward_rate_params[:,0] = np.log(forward_rate_params[:,0])
        if self.opt_type == 'params':
            x = x.reshape(len(self.reac_list), 3)
            for i, number in enumerate(self.reac_list):
                forward_rate_params[number] = np.multiply(forward_rate_params[number], x[i])
        #             forward_rate_params[number][1] = forward_rate_params[number][1] * x[3*i+1]
        #             forward_rate_params[number][2] = forward_rate_params[number][2] * x[3*i+2]

        sim_res = []
        temp_array = [self.temp]
        temp_array = np.array(temp_array).flatten()
        for ind, temp in enumerate(np.array(temp_array)):
            if self.opt_type == 'params':
                forward_rate_val = \
                    [np.exp(a[0]) * temp ** a[1] * np.exp((- a[2]) * self.conv / (GAS_CONST * temp))
                     for a in forward_rate_params]
            else:
                for i, parm in enumerate(self.reac_list):
                    forward_rate_val[parm] = \
                        np.multiply(np.power(10, x[i]), forward_rate_val[parm])

            if self.chemkin_data:
                chemkin = True
                free_energy_dict = generate_thermo_dict(self.energy,
                                                        self.smiles, temp)
                forward_rate_val = update_rate_constants_for_pressure(
                    self.chemkin_data, forward_rate_val, temp)
            else:
                chemkin = False
                free_energy_dict = build_free_energy_dict(self.energy, temp)

            rev_rate = build_reverse_rates(free_energy_dict,
                                           self.species_list, self.matrix,
                                           self.factor, forward_rate_val,
                                           temp, chemkin)

            time, sims = stiff_ode_solver(self.matrix, self.initial_y,
                                          forward_rate_val, rev_rate,
                                          self.third_body,
                                          sim_time=self.t_final,
                                          num_data_points=100)
            
            if time == 0:
                sims = self.opt_dist[ind] * 5.0
            sim_res.append(sims)

        if len(temp_array == 1):
            sim_res = sim_res[0]
        results, test_results = get_flatten_data(self.speciesnames, self.test_species, sim_res, self.temp,
                                                 self.sp_indices)
        simulations = [i / (0.1 * j) if j > 0.0 else i
                       for i, j in zip(results, self.opt_obs)]
        test_simulations = [i / (0.1 * j) if j > 0.0 else i for i, j in
                            zip(test_results, self.opt_test_obs)]
        #        print(simulations)
        return simulations, test_simulations

    def evaluation(self):
        #   observations = np.array(solve_ode_actual())

        observations = [j / (0.1 * j) if j > 0.0 else 0.0
                        for j in self.opt_obs]
        test_observations = [j / (0.1 * j) if j > 0.0 else 0.0
                             for j in self.opt_test_obs]
        return observations, test_observations

    def objectivefunction(self, simulation, evaluation):
        #         print(self.algorithm)
        if self.cost_func == 'sae':
            if self.algorithm in ['abc', 'fscabc']:
                objective_function = sae_func(evaluation[0], simulation[0])
                self.test_objectivefunction = sae_func(evaluation[1],
                                                       simulation[1])
            else:
                objective_function = - sae_func(evaluation[0], simulation[0])
                self.test_objectivefunction = - sae_func(evaluation[1],
                                                         simulation[1])
        else:
            if self.algorithm in ['abc', 'fscabc']:
                objective_function = getattr(spotpy.objectivefunctions,
                                             self.cost_func)(evaluation[0],
                                                             simulation[0])
                self.test_objectivefunction = \
                    getattr(spotpy.objectivefunctions,
                            self.cost_func)(evaluation[1], simulation[1])
            else:
                objective_function = - getattr(spotpy.objectivefunctions,
                                               self.cost_func)(evaluation[0],
                                                               simulation[0])
                self.test_objectivefunction = \
                    - getattr(spotpy.objectivefunctions,
                              self.cost_func)(evaluation[1], simulation[1])
        return objective_function

    def save(self, objectivefunctions, parameter, simulations, chains=None):
        parameter = list(parameter)
        line = str(objectivefunctions) + ',' + str(self.test_objectivefunction) \
               + ',' + str(parameter).strip('[]') + '\n'
        self.database.write(line)


def sae_func(predictions, targets):
    return ((np.array(predictions) - np.array(targets)) ** 2).sum()


def get_flatten_data(species, test_species, res, temp, sp_indeices):
    obs = []
    test_obs = []
    empty = 0
    temp_array = [temp]
    temp_array = np.array(temp_array).flatten()
    for i in range(len(temp_array)):
        if species == 'all':
            obs.append(np.array(res[i]).flatten())
            test_obs.append(empty)
        else:
            if len(temp_array == 1):
                dist = np.array(res)
            else:
                dist = np.array(res[i])
            for sp in species:
                obs.append(dist[:, sp_indeices[sp]])
            for sp in test_species:
                test_obs.append(dist[:, sp_indeices[sp]])
    opt_obs = np.array(obs).flatten()
    opt_test_obs = np.array(test_obs).flatten()

    return opt_obs, opt_test_obs


def create_parameter_name_list(reaction_id, par_type='rate'):
    if par_type == 'rate':
        total_list = ['k{}'.format(i) for i in reaction_id]
    else:
        a_list = ['A{}'.format(i) for i in reaction_id]
        n_list = ['n{}'.format(i) for i in reaction_id]
        e_list = ['E{}'.format(i) for i in reaction_id]
        # total_list = [a,  b, c for a_, b_, c_ in zip(a_list, n_list, e_list)]
        ini_list = [list(a) for a in zip(a_list, n_list, e_list)]
        total_list = [j for sub in ini_list for j in sub]

    return total_list


def optimization(pos, rep, reaction_list, opt_type, sp_indices, ini_val,
                 opt_dist, cost_function, forward_rate, rate_file,
                 energy_file, matrix, species_list, initial_y,
                 final_t, algorithm, temper, opt_species, third_body,
                 factor=1, chemkin_data=None, smiles=None, energy_conv='cal'):
    # print(sp_indices)
    parallel = "seq"
    dbformat = "custom"
    timeout = 1e6
    conver = 1
    if energy_conv == 'kcal':
        conver = KCAL_JL
    if energy_conv == 'cal':
        conver = CAL_JL
    if energy_conv == 'hartree':
        conver = HT_JL
    if energy_conv == 'KJ':
        conver = 1000
    if energy_conv == 'J':
        conver = 1
    spot_setup = ParamOptimize(reaction_list, opt_type, sp_indices, ini_val, opt_dist, cost_function,
                               forward_rate, rate_file, energy_file, matrix,
                               species_list, initial_y, final_t, third_body,
                               algorithm, temper, conver, opt_species,
                               factor=factor, chemkin_data=chemkin_data, smiles=smiles, pos=pos)

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


def multi_optimization(processes, rep, reac_list, opt_type, sp_indices,
                       ini_val, opt_dist, cost_function, forward_rate,
                       rate_file, energy_file, matrix, species_list,
                       initial_y, final_t, algorithms, temper,
                       opt_species='all', third_body=None,
                       chemkin_data=None, smiles=None, factor=1, energy_conv='cal'):
    print(processes)
    pool = mp.Pool(processes=processes)

    results = [pool.apply_async
               (optimization, args=(pos, rep, reac_list, opt_type, sp_indices,
                                    ini_val, opt_dist, cost_function,
                                    forward_rate, rate_file, energy_file,
                                    matrix, species_list, initial_y,
                                    final_t, al, temper, opt_species, third_body, factor,
                                    chemkin_data, smiles, energy_conv))
               for (pos, al) in enumerate(algorithms)]
    results = [p.get() for p in results]
    results.sort()  # to sort the results by input window width
    results = [r[1] for r in results]
    return results
