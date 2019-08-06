import numpy as np
from .constants import GAS_CONST, PR_ATM
from .constants import KCAL_JL, HT_JL
import math
import pandas as pd
import re
from itertools import islice


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


def set_paths_chemkin_files(my_path):
    """
    Set the absolute path to required files on the current machine.
    *** only required if using chemkin files****
    Parameters
    -------
    my_path                 : str
                            path where all the imput files are located
    Returns
    -------
    thermo_path            : str
                            path to the chemkin thermo file
    smile_path             : str
                            path to the file `species_smiles.dat`
    reactionlist_path      : str
                            path to the chemkin reaction file
    """

    thermo_path = my_path + '/thermo.dat'
    smile_path = my_path + '/species_smiles.dat'
    reactionlist_path = my_path + 'reaction.dat'
    return thermo_path, smile_path, reactionlist_path


def generate_smile_dict_chemkin(path):
    """
    creates a dictionary where key in the species formula
    and the value is the corresponding smiles
    *** only required if using chemkin files****
    paratemetrs
    -------
    path                   : str
                            path of the 'species_smiles.dat' file
    Returns
    -------
    dict_smile             : dict
                            dictionary key = formula, value = smiles

    """
    data = pd.read_csv(path)
    dict_smiles = dict(zip(data.Molecule, data.Smile))

    return dict_smiles


def parse_chemkin_thermo(file_name, dictionary):
    """
    parse the chemkin thermo_file, gets all the
    coefficient values for all the species
    Parameters
    -------
    file_name               : str
                            path of the 'thermo.dat' file
    dictionary              : dict
                            the smiles dictionary generated from
    Returns
    -------
    low                     : float
                            low temperature
    mid                     : float
                            mid temperature
    high                    : float
                            high temperature
    ***Each chemkin thermo file comes with a low, mid and high temperature
    values where the thermo coefficients are valid***
    dict_thermo_coeff       : dict
                            A dictionary where the key is the smiles of the
                            molecules and the values are the corresponding
                            coefficients from the chemkin thermo file
                            to calculate the thermodynamic properties
    """

    # number of lines in each chunk
    num = 4
    #     file = open(file_name, 'r')

    species_list = []
    thermo_coeff = []

    with open(file_name) as f:
        for i in range(2):
            line = f.readline()
            if i == 1:
                low = float(line.split()[0])
                mid = float(line.split()[1])
                high = float(line.split()[2])

        while True:
            next_n_lines = list(islice(f, num))
            if not next_n_lines:
                break
            smi, thermo = collect_thermo_params(next_n_lines, dictionary)
            species_list.append(smi)
            thermo_coeff.append(thermo)
    dict_thermo_coeff = dict(zip(species_list, thermo_coeff))
    #    print(thermo_dict)
    return low, mid, high, dict_thermo_coeff


def collect_thermo_params(lines, dictionary):
    """
    This fucntion generates a dictionary where key is the
    molecule smiles and values are lists of coefficients
    from the chemkin thermo file those are required to
    calculate the thermodynamic properties of each cpecies.
    Parameters
    -------
    lines               : list
                        chunk of four lines for each species
                        from chemkin thermo file
    dictionary          : dict
                        the smiles dictionary generated from
                        species_smiles.dat file
    Returns
    -------
    dictionary[species_name]   : str
                                smiles notation of the species
    thermo_list             : list
                            A list of coeffiecents from chemkin
                            thermo file
    """
    thermo_list = []
    species_name = None
    for i, line in enumerate(lines):
        if i == 0:
            species_name = line.split()[0]
        else:
            if i == 3:
                n = 60
            else:
                n = 75
            for start in range(0, n, 15):
                s = line[start:start + 15]
                thermo_list.append(float(s))
    #    print(thermo_params)
    #    print(dictionary[species_name])
    return dictionary[species_name], thermo_list


def generate_thermo_dict(file_path, smiles, temp):
    """
    This function generates the thrmodynamic properties
    of all the species present in the machanism from the
    chemkin thermo file.
    Parameters
    -------
    file_path           : str
                        path of the chemkin thermo file
    smiles              : dict
                        the smiles dictionary generated from
                        species_smiles.dat file
    temp                : float
                        system temperature
    Returns
    -------
    dic_thermo_values   : dict
                        a dictionary where the keys are the smiles of a
                        species and the values are a list of speciefic heat,
                        enthalpy, entropy and free energy of that epecies
                        at the given temperature
    """
    low_t, mid_t, high_t, thermo_params_dict = parse_chemkin_thermo(file_path, smiles)
    if temp < mid_t:
        i = 0
    else:
        i = 7

    species_names = []
    thermo_list = []
    for species, params in thermo_params_dict.items():
        #        print(species)
        species_names.append(species)
        specific_heat = (params[i + 0] + params[i + 1] * temp +
                         params[i + 2] * temp ** 2 +
                         params[i + 3] * temp ** 3 +
                         params[i + 4] * temp ** 4) * GAS_CONST
        enthalpy = (params[i + 0] + params[i + 1] * temp / 2 +
                    params[i + 2] * temp ** 2 / 3 +
                    params[i + 3] * temp ** 3 / 4 +
                    params[i + 4] * temp ** 4 / 5 +
                    params[i + 5] / temp) * temp * GAS_CONST
        entropy = (params[i + 0] * np.log(temp) + params[i + 1] * temp +
                   params[i + 2] * temp ** 2 / 2 +
                   params[i + 3] * temp ** 3 / 3 +
                   params[i + 4] * temp ** 4 / 4 + params[i + 6]) * GAS_CONST
        free_energy = enthalpy - entropy * temp
        list_entries = [specific_heat, enthalpy, entropy, free_energy]
        thermo_list.append(list_entries)
    #    print(species_names)
    dict_thermo_values = dict(zip(species_names, thermo_list))
    return dict_thermo_values


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
    return matrix, rate_final


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
                        conversion factor from given unit of energy to kJ
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

    return gibbs_energy_list, mol_change


def build_kmatrix_reverse(free_energy, species_list, stoic_mat,
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
                        conversion factor from given unit of energy to kJ
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

    gibbs_energy, change_mol = build_free_energy_change(free_energy,
                                                        species_list,
                                                        stoic_mat, factor,
                                                        chemkin)

    equilibrium_constants = [np.exp(-n * 1000/(GAS_CONST * temp))
                             for n in gibbs_energy]

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


def build_reaction_eqn(matrix, reaction_rate, rev_rate):
    """
    Builds the ode for single reactions
    matrix      : float
                the stoichiometric coefficient of the element
                for the specific reaction
    reaction_rate   : string
                    Reaction rate for the reaction calculated using
                    build_stoich_matrix function
    rev_rate         : 'yes' if reverse reactions are considered
                       'No' if reverse reactions are NOT considered
    Returns
    ----------
    rate        : string
                The rate equation for a specific reaction
                relative to a particular species

    """
    if matrix > 0:
        sign1 = ' + '
        sign2 = ' - '
    else:
        sign1 = ' - '
        sign2 = ' + '

    if abs(matrix) == 1.0:
        rate = sign1 + '%s' % re.split('-', reaction_rate)[0]
        if rev_rate == 'yes':
            rate += sign2 + '%s' % re.split('-', reaction_rate)[1]
        else:
            pass
    else:
        rate = sign1 + '%1.2f * %s' % (matrix, re.split('-', reaction_rate)[0])
        if rev_rate == 'yes':
            rate += sign2 + '%1.2f *%s' % (matrix, re.split('-', reaction_rate)[1])

        else:
            pass
    # print(rate)
    return rate


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
                                                species_list,
                                                species_indices)
    dydt_expressions = []
    for species in species_list:
        rate_equation = 'dy[%i]/dt = ' % (species_indices[species])
        i = species_indices[species]
        for j in range(len(rate_reac)):
            if abs(reac_matrix[j][i]) > 0:
                rate_equation += build_reaction_eqn(reac_matrix[j][i],
                                                    rate_reac[j], rev_rate)
            else:
                pass

        dydt_expressions.append(rate_equation)

    return dydt_expressions
