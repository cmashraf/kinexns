import numpy as np
from .constants import GAS_CONST, PR_ATM
from .constants import KCAL_JL, HT_JL, CAL_JL
import math
import pandas as pd
# import re
# from itertools import islice


def set_paths(my_path):

    """
    Set the absolute path to required files on the current machine.
    Parameters
    -------
    my_path                 : str
                            path where all the imput files are located
    Returns
    -------
    reactionlist_path     : str
                            path to the file `complete_reactionlist.dat`
    rateconstantlist_path : str
                            path to the file `complete_rateconstantlist.dat`
    free_energy_path    : str
                            path to the file 'free_energy_library.dat'
    """

    reactionlist_path = my_path + '/data/complete_reaction_list.dat'
    rateconstantlist_path = my_path + '/data/complete_rateconstant_list.dat'
    free_energy_path = my_path + '/data/free_energy_library.dat'

    return reactionlist_path, rateconstantlist_path, free_energy_path


class Reaction(object):

    """
    This is reaction class - it reads the reaction file and generate
    reactant list, product list and list of unique species in the mechanism.
    """

    def __init__(self):

        # initiating the class
        self.reactants_names = []
        self.products_names = []
        self.uniqueSpeciesList = []
        # species_names = []

    def get_reactants_name(self, line):
        """getting the reactants for each reaction
        Parameters
        ____________
        line        : str
                    line from the files
        Returs
        ____________
        reactants_names:    list
                            A list with reactant names and their
                            stoichiometric ratios in the reaction
        """

        for spec in line.split(','):
            if float(spec.split('_')[0].split()[0]) < 0:
                self.reactants_names.append((spec.split('_')[0].split()[0],
                                             spec.split('_')[1].split()[0]))
            # print(self.species_names)
        return self.reactants_names

    def get_products_name(self, line):
        """getting the products for each reaction
        Parameters
        ____________
        line        : str
                    line from the files
        Returs
        ____________
        reactants_names:    list
                            A list with product names and their
                            stoichiometric ratios in the reaction
        """

        for spec in line.split(','):
            if float(spec.split('_')[0].split()[0]) > 0:
                self.products_names.append((spec.split('_')[0].split()[0],
                                            spec.split('_')[1].split()[0]))
            # print(self.species_names)
        return self.products_names

    def unique_species_name(self, line, specieslist):
        """building the unique species list
        Parameters
        ____________
        line        : str
                    line from the files
        specieslist :     list
                            A list of species already in the mechanism
        Returs
        ____________
        reactants_names:    list
                            A list with reactant names and their
                            stoichiometric ratios in the reaction
        """

        # self.uniqueSpeciesList = species_list
        for spec in line.split(','):
            # self.uniqueSpeciesList = species_list
            #  If the species has already been added to the list then move on.
            if spec.split('_')[1].split()[0] in specieslist:
                self.uniqueSpeciesList = specieslist
                continue
            else:
                # print(self.uniqueSpeciesList)
                self.uniqueSpeciesList = specieslist
                self.uniqueSpeciesList.append(spec.split('_')[1].split()[0])
            # print(spec.split('_')[1].split()[0])
        return self.uniqueSpeciesList


def build_species_list(reaction_file):
    """
    Build reactnat and product list for each reaction. Also builds a list
    of unique species in the mechanism
    Parameters
    ----------
    reaction_file       : str
                           path to the file `complete_reaction_list.dat`
    Returns
    ----------
    reactant_list       : list
                         a list of the reactants and their stoichiometric
                         coeffs for each reaction
    product_list        : list
                         a list of the products and their stoichiometric
                         coeffs for each reaction
    species_list        : list
                        a list of unique species in the mechanism
    """

    # initializing reactant, product and unique species list
    reactant_list = []
    product_list = []
    species_name = []
    species_list = []
    file = open(reaction_file, 'r')

    for line in file:
        reac = Reaction()
        reactant_list.append(reac.get_reactants_name(line))
        product_list.append(reac.get_products_name(line))
        current_species = species_name
        # print(current_species)
        species_list = reac.unique_species_name(line, current_species)
        # print(species_name)

    species_list.sort()
    file.close()
    complete_list = [react + prod for react, prod in zip
                     (reactant_list, product_list)]
    species_indices = {species_list[i]: i for i in
                       range(0, len(species_list))}
    return complete_list, species_indices, species_list


def update_eff_dict(chemkin_data, species_list):
    """
    This function updates the third body efficiency dict generated
    by parsing chemkin reaction mechanism file. Since not all the
    species present in the efficienct dictionary might not be present
    in the mechanism, this function goes through the entire dictionary
    and deletes the species not present in mechanism
    Parameters
    ____________
    chemkin_data        :
                        the output generated from parse_chemkin_reaction
                        function
    species_list        : list
                        a list of unique species in the mechanism
    Returns
    ____________
   updated_eff_dictlist : list of dictionaries
                        a list of dictionaries with the species present in
                        the mechanism only
    """
    updated_eff_dictlist = []
    dictlist = chemkin_data[1]
    for i in range(len(dictlist)):
        keys = list(dictlist[i].keys())
        keys = [item for item in keys if item in species_list]
        values = [dictlist[i][k] for k in keys]
        updated_eff_dict = dict(zip(keys, values))
        updated_eff_dictlist.append(updated_eff_dict)
    return updated_eff_dictlist


def build_third_body_mat(chemkin_data, complete_list, species_list):
    """
    builds a 2D array like stoichiometric matrix where the number of
    row is equal to number of reactions and number of columns is equal
    to number of species. Each entry of the matrix corresponds the
    third body efficiency of that species in the reaction number
    corresponding to the row number
    Parameters
    ____________
    chemkin_data        : list of lists
                        the output generated from parse_chemkin_reaction
                        function
    complete_list       : list
                        A list of all the reactions with reactant and
                        product species and their stoichimetric coeffs
    species_list        : list
                        a list of unique species in the mechanism
    Retruns
    ____________
    third_body           : array
                         a 2D numpy array with third body efficiencies
    """

    eff_dict = update_eff_dict(chemkin_data, species_list)
    third_body = np.zeros((len(complete_list), len(species_list)), dtype=float)
    reaction_numbers = chemkin_data[0]
    el_numbers = []
    values = []
    for i, j in zip(range(len(eff_dict)), reaction_numbers):
        keys_list = list(eff_dict[i].keys())
        keys_list.sort()
        el_numbers.append([(j, species_list.index(k)) for k in keys_list])
        values.append(list(eff_dict[i].values()))
    #     #third_body_matrix[1][1]
    for row in reaction_numbers:
        third_body[row, :] = 1

    #    arr_el_numbers = np.array(el_numbers).reshape(-1, 2)
    arr_el_numbers = [x for sublist in el_numbers for x in sublist]
    arr_values = [x for sublist in values for x in sublist]

    for (i, j), val in zip(arr_el_numbers, arr_values):
        #        print(val)
        third_body[i, j] = val

    return third_body


def get_forward_rate_constants(parameters, temp, convert):
    """
    Generating the forward rate constants for each reaction
    Parameters:
    ____________
    parameters          : list
                        A list of Arrhenius paramters
    T                   : float, temperature
    convert             : str
                        unit conversion from 'convert' to JL
    Returns
    ____________
    forward_rate_parms  :    list
                        A list of forward rate constants (k_matrix)
    """
    factor = 0
    convert_val = {'cal': CAL_JL, 'kcal': KCAL_JL,
                   'hartrees': HT_JL, 'KJ': 1000, 'J': 1}
    factor = convert_val.get(energy_conv)

    forward_rates = (eval(parameters[0]) * temp ** eval(parameters[1]) *
                     np.exp((- eval(parameters[2]) * factor / (GAS_CONST * temp))))
    return forward_rates


class KineticParams(object):
    """
    This is the kinetic params class, they read the rates constant file,
    and generate the rate constants from the Arrhenius equations
    """

    def __init__(self):
        self.forward_rate_params = []
        self.forward_rates = []

    def get_forward_rate_parameters(self, line):
        """
        Reading the parameter file and parsing the useful infos
        Parameters
        ____________
        line        : str
                    line from the files
        Returns
        ____________
        forward_rate_parms  :    list
                            A list of Arrhenius paramters for the
                            forward reaction

        """

        self.forward_rate_params = [line.split(' ')[0], line.split(' ')[1],
                                    line.split(' ')[2].split()[0]]

        return self.forward_rate_params


def build_forward_rates(rateconstantlist, temp, convert='cal'):
    """
    This function builds the forward rate values for all
    the reactions. The Arrhenius rate parameters are found
    in the rateconstantlist file.
    Parameters
    ----------
    rateconstantlist        : str
                            path or name of the file to read from
    temp                    : float
                            temperature
    convert                 : str
                            unit conversion from 'convert' to JL
                            default = 'cal
    Returns
    ----------
    rate_constants      : list
                        a list of forward rate constants for all
                        the reactions
    """
    rate_constants = []
    params_list = []
    file = open(rateconstantlist, 'r')
    for line in file:
        f_params = KineticParams()
        params = f_params.get_forward_rate_parameters(line)
        params_list.append(params)
        rate_constants.append(get_forward_rate_constants(params, temp, convert))
    file.close()
    params_list = np.asarray(params_list)
    params_list = np.asfarray(params_list, float)
    return rate_constants, params_list


def update_rate_constants_for_pressure(chemkin_data, rate_constants, temp):
    """
    This function updates the forward rate coefficents of the pressure
    dependent reactions. It uses Troe formula to update the rate constants
    of the fall off reactions.
    Parameters
    ____________
    chemkin_data        : list of lists
                        the output generated from parse_chemkin_reaction
                        function
    rate_constants      : list
                        forward rate constants
    temp                    : float
                            temperature
    Returns
    ----------
    rate_constants      : list
                        a list of updated forward rate constants for all
                        the reactions
    """
    reaction_numbers = chemkin_data[2][0]
    rate_params = chemkin_data[2][1]
    troe_params = np.asarray(chemkin_data[2][2])

    for i, num in enumerate(reaction_numbers):
        k_0 = rate_constants[i]
        k_inf = rate_params[i][0] * temp ** rate_params[i][1] * \
            np.exp(- rate_params[i][2] * CAL_JL / (GAS_CONST * temp))
        p_r = k_0 / k_inf
        if len(troe_params[i]) == 1:
            troe_value = 1
        else:
            if troe_params[i][2] == 0:
                troe_params[i][2] = 1e-30
            f_cent = (1 - troe_params[i][0]) * \
                np.exp(- temp / troe_params[i][1]) + troe_params[i][0] *\
                np.exp(- temp / troe_params[i][2]) + \
                np.exp(- troe_params[i][3] / temp)
            c = -0.4 - 0.67 * np.log10(f_cent)
            n = 0.75 - 1.27 * np.log10(f_cent)
            f1 = (np.log10(p_r) + c) / (n - 0.14 * (np.log10(p_r) + c))
            log_f = np.log10(f_cent)/ (1 + f1 ** 2)
            troe_value = 10 ** log_f
        rate_constants[num] = k_inf * (p_r / (1 + p_r)) * troe_value
    return rate_constants


def build_free_energy_dict(free_energy_path, temp):
    """
    Build a dictionary of free energy at a given temperature for all
    the species present in the mechanism. It reads the file free_energy_path
    which is basically a library of gibbs free energy correction at
    different molecules at different temperatures.
    Parameters
    ----------
    free_energy_path     : str
                           path to the file `free_energy_library.dat`
    temp                 : float
                           temperature to calculate free energy
    Returns
    -------
    free_energy.    : dict
                    a dictionary where keys are unique species and values
                    are free energy of species at a given temperature
                    build from free_energy_library.dat
    """

    df = pd.read_csv(free_energy_path, sep='\t')

    if "{}K".format(temp) in df.columns:
        df["Free Energy @{}K".format(temp)] = df['electronic_energy'] +\
                                              df["{}K".format(temp)]
    else:
        temp_low = math.floor(temp / 100.0) * 100
        temp_high = math.ceil(temp / 100.0) * 100
        df["{}K".format(temp)] = ((df["{}K".format(temp_high)] -
                                  df["{}K".format(temp_low)]) *
                                  (temp - temp_low) / (temp_high - temp_low) +
                                  df["{}K".format(temp_low)])
        df["Free Energy @{}K".format(temp)] = (df['electronic_energy'] +
                                               df["{}K".format(temp)])
# print(df.head())

    free_energy = dict([(i, a) for i, a in
                        zip(df.smiles, df["Free Energy @{}K".format(temp)])])

    return free_energy


def build_stoic_matrix(complete_list, species_list, indices_species):
    """
    builds the stoichiometric matrix for the entire mechanism and
    then builds the rate equations for each reaction.
    Parameters
    ----------
    complete_list    : list
                      A list of all the reactions with reactant and
                      product species and their stoichimetric coeffs
    species_list     : list
                     A list of unique species in the mechanism
    indices_species  : dict
                     the dictionary species_indices
    Returns
    ----------
    matrix           : matrix
                     stoichiometric matrix for the entire
                     reaction mechanism
    rate_final       : str
                     list of strings of rate equations for each
                        reaction
    """

    matrix = np.zeros((len(complete_list), len(species_list)), dtype=float)
    rate_final = []
    for rxnindex, reac_list in enumerate(complete_list):
        rate = build_rate(reac_list, species_list, matrix, rxnindex, indices_species)
        rate_final.append(rate)
    return matrix


def build_free_energy_change(free_energy, species_list, stoic_mat, factor, chemkin=True):
    """
    Calculate the free energy changes for all the reactions
            delG = G(products) - G(reactanat)
    This is calculated from the complete lists of reactions
    and free_energy_dict
    Parameters
    ----------
    free_energy          : dict
                          A dictionary of free energies of all the species
                          at a given temperature, obtained from
                          build_free_energy_dict function or from chemkin
                          thermo files
    species_list        : list
                        a list of unique species in the mechanism
    stoic_mat           : nampy array
                        stoichiometric matrix of the mechanism
    factor              : float
                        conversion factor from given unit of energy to J
    chemkin             : bool
                        indicates if chemkin files are read as input files
                        default = True

    Returns
    -------
    gibbs_enenrgy_change : list
                         A list of free energy change for each reaction
    mol_change           : list
                         A list of (n_products - n_reactants)
                         for each reation
    """
    free_energy_sorted = {}
    #    gibs_energy_list = []
    for i in species_list:
        if chemkin:
            free_energy_sorted[i] = free_energy[i][3]
        else:
            free_energy_sorted[i] = free_energy[i]

    mol_change = stoic_mat.sum(axis=1, dtype=float)
    free_energy_list = list(free_energy_sorted.values())
    gibbs_energy_list = np.dot(stoic_mat, free_energy_list) * factor

    return gibbs_energy_list


def build_reverse_rates(free_energy, species_list, stoic_mat,
                        factor, forward_rates, temp, chemkin=True):
    """"
    Calculates the reverse rate constants for all the reactions
    using the free energy change through the following steps
    1. Use delG from build_free_energy_change
    to calculate the equlilibrium constant
    Keq = exp (- delG/Gas Const * temp)
    2. Use the following equation to calculate the reverse rate constant
    Keq = Kf / Kr * (Gas Const * temp / Pressure)^n
    where n = total number of product molecules -
    total number of reactant molecules
    Parameters
    ----------
    free_energy          : dict
                          A dictionary of free energies of all the species
                          at a given temperature, obtained from
                          build_free_energy_dict function or from chemkin
                          thermo files
    species_list        : list
                        a list of unique species in the mechanism
    stoic_mat           : nampy array
                        stoichiometric matrix of the mechanism
    factor              : float
                        conversion factor from given unit of energy to J
    chemkin             : bool
                        indicates if chemkin files are read as input files
                        default = True
    forward_rates        : A list of forward rate constants for all the
                         reactions obtained from build_forward_reaction_rates
    temp                : float
                           temperature to calculate free energy
    Returns
    -------
    reverse_rates       : list
                         A list of reverse rate constants
    """
    gibbs_energy = build_free_energy_change(free_energy, species_list,
                                            stoic_mat, factor, chemkin=chemkin)
    change_mol = stoic_mat.sum(axis=1, dtype=float)
    equilibrium_constants = [np.exp(-n / (GAS_CONST * temp))
                             for n in gibbs_energy]
    #    print(forward_rates)
    reverse_rates = [(a / b) * (GAS_CONST * temp * 1000 / PR_ATM) ** c
                     if c < 3 else 0 for (a, b, c) in
                     zip(forward_rates, equilibrium_constants, change_mol)]

    return reverse_rates


def build_concentartion(species_id, number):
    """
    Builds the concentration component for each species
    Parameters
    ----------
    species_id      : int
                     species id from the species_indices dictionary
    number          : float
                    stoichiometric co-eff of the species
                    for specific reactions
    Returns
    ----------
    concentration   : string
                    the concentration component of that particular
                    species
    """

    if abs(float(number)) == 1:
        concentration = '* y[%s] ' % species_id
    else:
        concentration = '* y[%s] ** %s ' % (species_id, abs(float(number)))

    return concentration


def build_rate(reac_prod, spc_list, matrix, index, indices_species):
    """
    Build the rate equation for each species
    Parameters
    ----------
    reac_prod       : list
                    list of the reactants and products
                    with stoiciometric coeffs for each reaction
    spc_list        : list
                    A list of unique species in the mechanism
    matrix          : matrix
                    the stoichiometric matrix that is being built
                    in build_stoic_matri function
    index           : int
                    index of the reaction
    indices_species  : dict
                     the dictionary species_indices
    Returns
    ----------
    rate_reac       : str
                    rate equation for individual reactions
    """
    concentration_f = ''
    concentration_r = ''
    rate_f = 'kf[%s] ' % index
    rate_r = '- kr[%s]' % index
    for x in range(len(reac_prod)):
        species = reac_prod[x][1]
        for i in range(len(spc_list)):
            if i == indices_species[species]:
                matrix[index][i] = float(reac_prod[x][0])
                if float(reac_prod[x][0]) < 0:
                    concentration_f += build_concentartion(i, reac_prod[x][0])
                else:
                    concentration_r += build_concentartion(i, reac_prod[x][0])

    rate_reac = rate_f + concentration_f + rate_r + concentration_r

    return rate_reac
#
#
# def build_reaction_eqn(matrix, reaction_rate, rev_rate):
#     """
#     Builds the ode for single reactions
#     matrix      : float
#                 the stoichiometric coefficient of the element
#                 for the specific reaction
#     reaction_rate   : string
#                     Reaction rate for the reaction calculated using
#                     build_stoich_matrix function
#     rev_rate         : 'yes' if reverse reactions are considered
#                        'No' if reverse reactions are NOT considered
#     Returns
#     ----------
#     rate        : string
#                 The rate equation for a specific reaction
#                 relative to a particular species
#
#     """
#     if matrix > 0:
#         sign1 = ' + '
#         sign2 = ' - '
#     else:
#         sign1 = ' - '
#         sign2 = ' + '
#
#     if abs(matrix) == 1.0:
#         rate = sign1 + '%s' % re.split('-', reaction_rate)[0]
#         if rev_rate == 'yes':
#             rate += sign2 + '%s' % re.split('-', reaction_rate)[1]
#         else:
#             pass
#     else:
#         rate = sign1 + '%1.2f * %s' % (matrix, re.split('-', reaction_rate)[0])
#         if rev_rate == 'yes':
#             rate += sign2 + '%1.2f *%s' % (matrix, re.split('-', reaction_rate)[1])
#
#         else:
#             pass
#     # print(rate)
#     return rate
#
#
# def build_dydt_list(complete_list, species_list,
#                     species_indices, rev_rate='yes'):
#     """
#     builds the stoichiometric matrix for the entire mechanism and
#     then builds the rate equations for each reaction.
#     Parameters
#     ----------
#     complete_list    : list
#                       A list of all the reactions with reactant and
#                       product species and their stoichimetric coeffs
#     species_list     : list
#                      A list of unique species in the mechanism
#     species_indices  : dict
#                      the dictionary speciesindices
#     rev_rate         : 'yes' if reverse reactions are considered
#                        'No' if reverse reactions are NOT considered
#     Returns
#     ----------
#     dydt_expressions : list
#                      A lsit of all the differential equations
#                      related to the mechanism
#     """
#
#     reac_matrix, rate_reac = build_stoic_matrix(complete_list,
#                                                 species_list,
#                                                 species_indices)
#     dydt_expressions = []
#     for species in species_list:
#         rate_equation = 'dy[%i]/dt = ' % (species_indices[species])
#         i = species_indices[species]
#         for j in range(len(rate_reac)):
#             if abs(reac_matrix[j][i]) > 0:
#                 rate_equation += build_reaction_eqn(reac_matrix[j][i],
#                                                     rate_reac[j], rev_rate)
#             else:
#                 pass
#
#         dydt_expressions.append(rate_equation)
#
#     return dydt_expressions
