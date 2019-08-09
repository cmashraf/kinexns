import numpy as np
from .constants import GAS_CONST, PR_ATM
from .constants import KCAL_JL, HT_JL, CAL_JL
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
    reactionlist_path = my_path + '/reaction.dat'
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


def species_string(string_lists, smiles, reactant=True):
    """
    This function generates the reactant and the product strings that
    to write in the reaction file from the chemkin reaction file to make
    the reaction strings compatible for this package.
    Parameters
    -------
    string_lists        : list of strings
                        the reactants and products list parsed from
                        chemkin reaction file
    smiles              : dict
                        smiles dictionary of all the molecules
    reactant            : bool
                        indictaes if the list is a reactant or product list
                        default = True indicates reactant list
    Returns
    -------
    formatted_string    : str
                        species smiles and stoichiometric coeff formatted
                        in proper way to be written in the file
                        format: coeff_smiles
    """

    formatted_string = []
    for i in range(len(string_lists)):
        #        print(s_prod)

        if ' ' in string_lists[i]:
            stoic = string_lists[i].split()[0]
            species = smiles[string_lists[i].split()[1]]
            # print(species)
        else:
            stoic = '1.0'
            species = smiles[string_lists[i]]
            # print(species)

        if reactant:
            stoic_val = -1 * float(stoic)
            stoic = str(stoic_val)
        formatted_string.append(stoic + '_' + species)

    formatted_string = ','.join(formatted_string)
    return formatted_string


def write_reactions(string, smiles, file_reaction, file_rate_cons):
    """
    This function gets a reaction line from chemkin reaction file,
    and then formats it to be written in a format that could be recognised
    by the package.
    Parameters
    -------
    string          : str
                    line that has been identified as a reaction line from
                    chemkin reaction file.
    smiles          : dict
                    smiles dictionary of all the molecules
    file_reaction   : textIO
                    opened file name or path where the reactions
                    will be written
    file_rate_cons       :textIO
                    opened file name or path where the reaction rates
                    will be written
    Returns
    -------
    """
    s = string.split('<=>')
    s2 = s[1].split()
    s_reac = s[0].split('+')
    s_pr = ' '.join(s2[:-3])
    #            print(s_pr)
    s_prod = s_pr.split('+')
    #            print(s_prod)
    reac_string = species_string(s_reac, smiles)
    prod_string = species_string(s_prod, smiles, reactant=False)

    reaction_string = reac_string + ',' + prod_string
    #    print(reaction_string)
    file_reaction.write(reaction_string + '\n')
    file_rate_cons.write(s2[-3] + ' ' + s2[-2] + ' ' + s2[-1] + '\n')


def parse_chemkin_reaction(file_name, smiles, file_reac, file_rate):
    """
    This function parse the chemkin reaction mechanism file with the
    help of two other functions (write_reactions and species_string and
    writes two seperate files with reactions and associated rate constants
    that will be recognised by this package.
    Parameters
    -------
    file_name       : str
                    name or path of the chemkin reaction mechanism file
    smiles          : dict
                    smiles dictionary of all the molecules
    file_reac       : str
                    file name or path where the reactions will be written
    file_rate       : str
                    file name or path where the reaction rates
                    will be written
    Returns
    -------
    three_body_reactions    : list
                            A list of rection numbers those have third body effects
    three_body_eff          : list of dictionaries
                            each list contains a dictionary with simes as
                            keys and third body efficiencies as values
                            for the species those have a tbe other than
                            one corresponding to the reaction numbers in
                            third_body_reactions
    pr_dependent            : lists of lists
                            three seperate lists
                            first list is a list of reaction numbers
                            those have pressure depedence
                            second list is a list of lists where each
                            sublist contains Arrhenius parameters to
                            calculate k0 values for the reaction
                            number in the first list
                            third list is a list of lists where each
                            sublist contains Troe parameters to calculate
                            troe function values for the
                            reaction number in the first list
    """

    file = open(file_name, 'r')
    file_reactions = open(file_reac, 'w+')
    file_rate_constants = open(file_rate, 'w+')
    i = 0
    three_body_reactions = []
    pr_dependent_reactions = []
    k_low = []
    troe_value = []
    three_body_eff = []
    low_count = 0
    troe_count = 0
    for line in file:
        #        print(line)

        if '+M' in line:
            line = re.sub(r'\+M', '', line)
            three_body_reactions.append(i)

        if '()' in line:
            line = re.sub(r'[()]', '', line)
            pr_dependent_reactions.append(i)
        if any(item in line for item in ['/', 'LOW', 'TROE']):
            if 'LOW' in line:
                low_count = low_count + 1
                s = line.split()
                k_low.append([float(s[2]), float(s[3]), float(s[4])])
            elif 'TROE' in line:
                troe_count = troe_count + 1
                if troe_count != low_count:
                    troe_value.append([0])
                    troe_count = troe_count + 1
                s = line.split()
                if len(s) == 6:
                    troe_value.append([float(s[1]), float(s[2]),
                                       float(s[3]), float(s[4])])
                else:
                    troe_value.append([float(s[1]), float(s[2]),
                                       float(s[3]), 1e10])
            else:
                s = line.split()
                keys = [s[i].split('/')[0] for i in range(len(s))]
                keys_smi = [smiles[key] for key in keys]
                values = [s[i].split('/')[1] for i in range(len(s))]

                dictionary = dict(zip(keys_smi, values))
                #                dictList.append()
                three_body_eff.append(dictionary)

        if '=' in line:
            write_reactions(line, smiles, file_reactions, file_rate_constants)
            i = i + 1

    pr_dependent = [pr_dependent_reactions, k_low, troe_value]

    file_reactions.close()
    file_rate_constants.close()

    #    print(three_body_eff)
    return three_body_reactions, three_body_eff, pr_dependent


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

    def get_forward_rate_constants(self, parameters, temp, convert):
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
        if convert == 'cal':
            factor = CAL_JL
        if convert == 'kcal':
            factor = KCAL_JL
        self.forward_rates = (eval(parameters[0]) * temp ** eval(parameters[1]) *
                              np.exp((- eval(parameters[2]) * factor / (GAS_CONST * temp))))
        return self.forward_rates


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
    for line in open(rateconstantlist, 'r').readlines():
        f_params = KineticParams()
        params = f_params.get_forward_rate_parameters(line)
        rate_constants.append(f_params.get_forward_rate_constants(params, temp, convert))

    return rate_constants


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
            troe_value = (1 - troe_params[i][0]) * \
                np.exp(- temp / troe_params[i][1]) + troe_params[i][0] *\
                np.exp(- temp / troe_params[i][2]) + \
                np.exp(- troe_params[i][3] / temp)
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
    gibbs_energy = build_free_energy_change(free_energy, species_list, stoic_mat, 0.001, chemkin=True)
    change_mol = stoic_mat.sum(axis=1, dtype=float)
    equilibrium_constants = [np.exp(-n * 1000 / (GAS_CONST * temp))
                             for n in gibbs_energy]
    #    print(forward_rates)
    reverse_rates = [(a / b) * (GAS_CONST * temp * 1000 / PR_ATM) ** c
                     for (a, b, c) in
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
