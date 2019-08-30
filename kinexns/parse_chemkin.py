import numpy as np
from .constants import GAS_CONST, PR_ATM
from .constants import KCAL_JL, HT_JL, CAL_JL
import math
import pandas as pd
import re
from itertools import islice


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
    low, mid, high = 0.0, 0.0, 0.0
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


def parse_pr_three_body_data(line_data, k_low, troe_value, three_body_eff, low_count, troe_count, smiles):
    if any(item in line_data for item in ['/', 'LOW', 'TROE']):
        if 'LOW' in line_data:
            low_count = low_count + 1
            s = line_data.split()
            k_low.append([float(s[2]), float(s[3]), float(s[4])])
        elif 'TROE' in line_data:
            troe_count = troe_count + 1
            if troe_count != low_count:
                troe_value.append([0])
                troe_count = troe_count + 1
            s = line_data.split()
            if len(s) == 6:
                troe_value.append([float(s[1]), float(s[2]), float(s[3]), float(s[4])])
            else:
                troe_value.append([float(s[1]), float(s[2]), float(s[3]), 1e10])
        else:
            s = line_data.split()
            keys = [s[i].split('/')[0] for i in range(len(s))]
            keys_smi = [smiles[key] for key in keys]
            values = [s[i].split('/')[1] for i in range(len(s))]

            dictionary = dict(zip(keys_smi, values))
            #                dictList.append()
            three_body_eff.append(dictionary)
    return k_low, troe_value, three_body_eff, low_count, troe_count


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

        if '=' in line:
            write_reactions(line, smiles, file_reactions, file_rate_constants)
            i = i + 1
        k_low, troe_value, three_body_eff, low_count, troe_count = \
            parse_pr_three_body_data(line, k_low, troe_value, three_body_eff,
                                     low_count, troe_count, smiles)
    pr_dependent = [pr_dependent_reactions, k_low, troe_value]

    file.close()
    file_reactions.close()
    file_rate_constants.close()

    #    print(three_body_eff)
    return three_body_reactions, three_body_eff, pr_dependent
