import numpy as np
from .constants import GAS_CONST, PR_ATM
from .constants import KCAL_JL, HT_JL
import math
import pandas as pd
import re


def set_paths(my_path):

    """
    Set the absolute path to required files on the current machine.

    Returns
    -------
    reactionlist_path     : str
                            path to the file `complete_reactionlist.dat`
    rateconstantlist_path : str
                            path to the file `complete_rateconstantlist.dat`
    compositionlist_path  : str
                            path to the file `compositionlist.dat`
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

    for line in open(reaction_file, 'r').readlines():
        reac = Reaction()
        reactant_list.append(reac.get_reactants_name(line))
        product_list.append(reac.get_products_name(line))
        current_species = species_name
        # print(current_species)
        species_list = reac.unique_species_name(line, current_species)
        # print(species_name)

    species_list.sort()

    return reactant_list, product_list, species_list


class KineticParams(object):
    """
    This is the kinetic params class, they read the rates constant file,
    and generate the rate constants from the Arrhenius equations
    """

    def __init__(self):
        self.forward_rate_params = []
        self.forward_rates = []
        # self.forward_E = []
        # self.uniqueSpeciesList = []
        # species_names = []

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

    def get_forward_rate_constants(self, parameters, temp):
        """
        Generating the forward rate constants for each reaction
        Parameters
        ____________
        parameters          : list
                            A list of Arrhenius paramters
        T                   : float, temperature
        Returns
        ____________
        forward_rate_parms  :    list
                            A list of forward rate constants (k_matrix)
        """

        self.forward_rates = (eval(parameters[0]) *
                              np.exp(- eval(parameters[2]) *
                              KCAL_JL / (GAS_CONST * temp)))
        return self.forward_rates


def build_kmatrix_forward(rateconstantlist, temp):

    rate_constants = []
    for line in open(rateconstantlist, 'r').readlines():
        f_params = KineticParams()
        params = f_params.get_forward_rate_parameters(line)
        rate_constants.append(f_params.get_forward_rate_constants(params, temp))

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
        df["{}K".format(temp)] = (df["{}K".format(temp_high)] -
                                  df["{}K".format(temp_low)]) *\
                                 (temp - temp_low) / (temp_high - temp_low)\
                                 + df["{}K".format(temp_low)]
        df["Free Energy @{}K".format(temp)] = df['electronic_energy'] + \
                                              df["{}K".format(temp)]
# print(df.head())

    free_energy = dict([(i, a) for i, a in
                        zip(df.smiles, df["Free Energy @{}K".format(temp)])])

    return free_energy


def build_free_energy_change(complete_list, free_energy):
    """
    Calculate the free energy changes for all the reactions
            delG = G(products) - G(reactanat)
    This is calculated from the complete lists of reactions
    and free_energy_dict
    Parameters
    ----------
    complete_list        : list
                           A list of all the reactions with reactant and
                         product species and their stoichimetric coeffs
    free_energy          : dict
                          A dictionary of free energies of all the species
                          at a given temperature, obtained from
                          build_free_energy_dict function
    Returns
    -------
    gibbs_enenrgy_change : list
                         A list of free energy change for each reaction
    mol_change           : list
                         A list of (n_products - n_reactants)
                         for each reation
    """

    mol_change = []
    gibbs_energy_list = []

    for i, item in enumerate(complete_list):
        n_reac = 0
        n_prod = 0
        reac_free_energy = 0
        prod_free_energy = 0
        for entry in item:

            if float(entry[0]) < 0:
                n_reac = n_reac + abs(float(entry[0]))
                reac_free_energy = (abs(float(entry[0])) *
                                    free_energy[entry[1]] +
                                    reac_free_energy)
            else:
                prod_free_energy = (abs(float(entry[0])) *
                                    free_energy[entry[1]]
                                    + prod_free_energy)
                n_prod = n_prod + abs(float(entry[0]))
        # print(n_reac)
        mol_change.append(n_prod - n_reac)
        # print(mol_change)
        gibbs_energy_list.append((prod_free_energy - reac_free_energy)
                                 * HT_JL)

    # equilibrium_constants = [np.exp(-n * 1000/(GAS_CONST * T))
    #                         for n in gibbs_energy_list]
    # reverse_rates = [(a / b) * 1000 * (GAS_CONST * T / PR_ATM) ** c
    #                 if c < 3 else 0 for (a, b, c) in zip(forward_rates,
    #                 equilibrium_constants, mol_change)]

    return gibbs_energy_list, mol_change


def build_kmatrix_reverse(complete_list, free_energy,
                          forward_rates, temp):
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
    complete_list        : list
                           A list of all the reactions with reactant and
                         product species and their stoichimetric coeffs
    free_energy          : dict
                          A dictionary of free energies of all the species
                          at a given temperature, obtained from
                          build_free_energy_dict function
    forward_rates        : A list of forward rate constants for all the
                         reactions obtained from build_forward_reaction_rates
    temp                : float
                           temperature to calculate free energy
    Returns
    -------
    reverse_rates       : list
                         A list of reverse rate constants
    """

    gibbs_energy, change_mol = build_free_energy_change(complete_list,
                                                        free_energy)

    equilibrium_constants = [np.exp(-n * 1000/(GAS_CONST * temp))
                             for n in gibbs_energy]

    reverse_rates = [(a / b) * (GAS_CONST * temp * 1000 / PR_ATM) ** c
                     if c < 3 else 0 for (a, b, c) in
                     zip(forward_rates, equilibrium_constants, change_mol)]

    return reverse_rates


def build_stoic_matrix(complete_list, species_list, species_indices):
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
    species_indices  : dict
                     the dictionary speciesindices
    Returns
    ----------
    matrix           : matrix
                     stoichiometric matrix for the entire
                     reaction mechanism
    rate_final       : str
                     list of strings of rate equations for each
                     reaction (both forward and reverse)
    """
    matrix = np.zeros((len(complete_list), len(species_list)), dtype=float)
    rate_final = []
    for rxnindex, reac_list in enumerate(complete_list):
        # rate = ''
        rate_f = 'kf[%s] ' % rxnindex
        rate_r = '- kr[%s]' % rxnindex
        concentration_f = ''
        concentration_r = ''
        for x in range(len(reac_list)):
            species = reac_list[x][1]
            for i in range(len(species_list)):
                if i == species_indices[species]:
                    matrix[rxnindex][i] = float(reac_list[x][0])
                    if float(reac_list[x][0]) < 0:
                        if abs(float(reac_list[x][0])) == 1:
                            concentration_f += '* y[%s] ' % i
                        else:
                            concentration_f += '* y[%s] ** %s ' \
                                               % (i, abs(float
                                                         (reac_list[x][0])))
                    else:
                        if abs(float(reac_list[x][0])) == 1:
                            concentration_r += '* y[%s] ' % i
                        else:
                            concentration_r += '* y[%s] ** %s ' \
                                               % (i, float(reac_list[x][0]))

        rate = rate_f + concentration_f + rate_r + concentration_r
        # print(rate)
        rate_final.append(rate)
    return matrix, rate_final


def build_dydt_list(complete_list, species_list,
                    species_indices, rev_rate='yes'):
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
    species_indices  : dict
                     the dictionary speciesindices
    rev_rate         : 'yes' if reverse reactions are considered
                       'No' if reverse reactions are NOT considered
    Returns
    ----------
    dydt_expressions : list
                     A lsit of all the differential equations
                     related to the mechanism
    """

    reac_matrix, rate_reac = build_stoic_matrix(complete_list,
                                                species_list, species_indices)
    dydt_expressions = []
    for species in species_list:
        rate_equation = 'dy[%i]/dt = ' % (species_indices[species])
        i = species_indices[species]
        for j in range(len(rate_reac)):
            if reac_matrix[j][i] > 0:
                if reac_matrix[j][i] == 1.0:
                    if rev_rate == 'yes':
                        rate_equation += '+ %s' % (rate_reac[j])
                    else:
                        rate_equation += '+ %s' % (re.split('-', rate_reac[j])[0])
                else:
                    rate_equation += '+ %1.2f * %s' % \
                                     (reac_matrix[j][i],
                                      re.split('-', rate_reac[j])[0])
                    if rev_rate == 'yes':
                        rate_equation += '- %1.2f * %s' % \
                                         (reac_matrix[j][i],
                                          re.split('-', rate_reac[j])[1])
                    else:
                        pass
            elif reac_matrix[j][i] < 0:
                if reac_matrix[j][i] == -1.:
                    rate_equation += '- %s' % (re.split('-', rate_reac[j])[0])
                    if rev_rate == 'yes':
                        rate_equation += '+ %s' % \
                                         (re.split('-', rate_reac[j])[1])
                    else:
                        pass
                else:
                    rate_equation += '- %1.2f * %s' % \
                                     (abs(reac_matrix[j][i]),
                                      re.split('-', rate_reac[j])[0])
                    if rev_rate == 'yes':
                        rate_equation += '+ %1.2f * %s' % \
                                         (abs(reac_matrix[j][i]),
                                          re.split('-', rate_reac[j])[1])
                    else:
                        pass
            else:
                pass
        dydt_expressions.append(rate_equation)

    return dydt_expressions
