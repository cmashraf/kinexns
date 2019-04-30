import os
import numpy as np
from .constants import GAS_CONST


def set_paths():
    
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
    
    #module_dir = os.getcwd().split('ligpy_utils')[0]
    myPath = os.path.dirname(os.path.abspath(__file__))
    reactionlist_path = myPath + '/data/complete_reaction_list.dat'
    rateconstantlist_path = myPath + '/data/complete_rateconstant_list.dat'
    compositionlist_path = myPath + '/data/compositionlist.dat'

    return reactionlist_path, rateconstantlist_path, compositionlist_path


class Reaction(object):
    
    """
    This is reaction class - it reads the reaction file and generate  
    reactant list, product list and list of unique species in the mechanism.
    """    
    
    def __init__(self):

        #initiating the class
        self.reactants_names = []
        self.products_names = []
        self.uniqueSpeciesList = []
        #species_names = []
        
   
    def getReactantsName(self, line):
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
            #print(self.species_names)
        return self.reactants_names
    
    def getProductsName(self, line):        
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
            #print(self.species_names)
        return self.products_names
    
    def uniqueSpeciesName(self, line, species_list):
        """building the unique species list
        Parameters
        ____________
        line        : str
                    line from the files
        species_list :     list
                            A list of species already in the mechanism
        Returs
        ____________
        reactants_names:    list
                            A list with reactant names and their 
                            stoichiometric ratios in the reaction
        """

        #self.uniqueSpeciesList = species_list
        for spec in line.split(','):
            #self.uniqueSpeciesList = species_list
            # If the species has already been added to the list then move on.
            if spec.split('_')[1].split()[0] in species_list:
                self.uniqueSpeciesList = species_list
                continue
            else:
                #print(self.uniqueSpeciesList)
                self.uniqueSpeciesList = species_list
                self.uniqueSpeciesList.append(spec.split('_')[1].split()[0])
            #print(spec.split('_')[1].split()[0])
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
    __________

    reactant_list       : list
                         a list of the reactants and their stoichiometric
                         coeffs for each reaction
    product_list        : list
                         a list of the products and their stoichiometric
                         coeffs for each reaction
    species_list        : list
                        a list of unique species in the mechanism
    """ 

    #initializing reactant, product and unique species list
    reactant_list = []
    product_list = []
    species_name = []

    for line in open(reaction_file, 'r').readlines():
        reac = Reaction()
        reactant_list.append(reac.getReactantsName(line))
        product_list.append(reac.getProductsName(line))
        current_species = species_name
        #print(current_species)
        species_list = reac.uniqueSpeciesName(line, current_species)
        #print(species_name)
    species_list.sort()

    return reactant_list, product_list, species_list



class Kinetic_params(object):
    """
    This is the kinetic params class, they read the rates constant file,
    and generate the rate constants from the Arrhenius equations
    """  

    def __init__(self):
        self.forward_rate_params = []
        self.forward_rates = []
        #self.forward_E = []
        #self.uniqueSpeciesList = []
        #species_names = []
    
    def getForwardRateParameters(self, line):
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
    
    def getForwardRateConstant(self, parameters, T):
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
        

        self.forward_rates = eval(parameters[0]) * np.exp(- eval(parameters[2])/
                                                               (GAS_CONST * T))
        return self.forward_rates


def build_kmatrix_forward(rateconstantlist, temp):
    
    rate_constants = []
    for line in open(rateconstantlist, 'r').readlines():
        f_params = Kinetic_params()
        params = f_params.getForwardRateParameters(line)
        rate_constants.append(f_params.getForwardRateConstant(params, temp))
    
    return rate_constants


def build_reac_prod_dict(reac_list, prod_list, speciesindices):
    """
    Build a dictionary of the reactants involved in each reaction,
    along with their stoichiometric coefficients.  The keys of the
    dictionary are the reaction numbers, the values are lists of lists
    [[reactant1index, -1*coeff1],...]
    Parameters
    ----------
    completereactionlist : str
                           path to the file `complete_reaction_list.dat`
    speciesindices       : dict
                           the dictionary speciesindices from
                           get_speciesindices()
    Returns
    -------
    reactant_dict : dict
                    a dictionary where keys are reaction numbers and values
                    are lists of lists with the reactant species id and their
                    stoichiometric coefficients for each reaction
    product_dict : dict
                    a dictionary where keys are reaction numbers and values
                    are lists of lists with the product species id and their
                    stoichiometric coefficients for each reaction
    """
    reactant_dict = {}
    for rxnindex, reaction in enumerate(reac_list):
        reactants = []
        #
        for x in range(len(reaction)):
            # if the species is a reactant
         #   if float(x.split('_')[0]) < 0:
            reactants.append([speciesindices[reaction[x][1]],
                                -1*float(reaction[x][0])])
            #    in preceding line: *-1 because I want the |stoich coeff|
        reactant_dict[rxnindex] = reactants
        
    products_dict = {}
    for rxnindex, reaction in enumerate(prod_list):
        products = []
        #
        for x in range(len(reaction)):
            # if the species is a reactant
         #   if float(x.split('_')[0]) < 0:
            products.append([speciesindices[reaction[x][1]],
                                1*float(reaction[x][0])])
            #    in preceding line: *-1 because I want the |stoich coeff|
        products_dict[rxnindex] = products
    return reactant_dict, products_dict


def build_reac_species_dict(reacprodlist, specieslist):
    """
    Build a dictionary where keys are species and values are lists with the
    reactions that species is involved in, that reaction's sign in the net
    rate equation, and the stoichiometric coefficient of the species in that
    reaction.
    Parameters
    ----------
    reacprodlist : list
                        a list of both reactants and products and their 
                        stoichiometric co-effs
    specieslist  : list
                        a list of unique species in the mecahnism
    
    Returns
    -------
    reac_species : dict
                   keys are the species in the model; values are lists of
                   [reaction that species is involved in,
                   sign of that species in the net rate equation,
                   sign for forward reaction ('-1', if reactant, '+1', if product),
                   sign for backward reaction ('+1', if reactant, '-1', if product)]
    """
    #specieslist = get_specieslist(set_paths()[0])
    reac_species = {}
    for species in specieslist:
        # This loop makes a list of which reactions "species" takes part in
        # and what sign that term in the net rate eqn has
        # and what the stoichiometric coefficient is
    
        reactions_involved = []
        for rxnindex, reac_list in enumerate (reacprodlist):
            for x in range(len(reac_list)):
                # If the species being iterated over is part of this reaction
                if species == reac_list[x][1]:
                    # if the species is a reactant
                    if float(reac_list[x][0]) < 0:
                        reactions_involved.append(
                            [rxnindex, -1, str(-1), '+'+str(1)])
                    
                    # if the species is a product
                    if float(reac_list[x][0]) > 0:
                        reactions_involved.append(
                            [rxnindex, 1, '+'+str(1), str(-1)])
    
        reac_species[species] = reactions_involved
    return reac_species


def build_rate_eqn(k_mat, r_dict, human, forward):
    
    """ This function writes the list of rate expressions for each reaction.
    Parameters
    ----------
    kmat               : list
                         A list of reaction rate contstants (k_forward or k_reverse)
    r_dict             : dictionary
                         reactant or product directory
    human              : str, optional
                         indicate whether the output of this function should
                         be formatted for a human to read ('yes'). Default
                         is 'no'
    forward             : str
                        reaction type,if 'yes', it is forward reaction
                        default is 'yes'
    Returns
    -------
    rates_list : list
                A list of the rate expressions for all the reactions in the mecahnism
    """

    rates_list = []
    for i, line in enumerate(k_mat):
        if forward == 'yes':
            rate = 'rate_f[%s] = kf(T,%s) ' % (i, i)
        else:
            rate = 'rate_r[%s] = kf(T,%s) ' % (i, i)
        concentrations = ''
        for entry in r_dict[i]:
            if entry == 'n':   # if there is no reaction
                concentrations = '* 0'
                break
            else:
                if human == 'no':
                    concentrations += '* y[%s]**%s ' % (entry[0], entry[1])
                elif human == 'yes':
                    concentrations += '* [%s]**%s ' % \
                        (indices_to_species[entry[0]], entry[1])
                else:
                    raise ValueError('human must be a string: yes or no')
    
        rate += concentrations
        
        #rate = rate_reactant + rate_product
        rates_list.append(rate)
        
    return rates_list
    


def build_rates_list(reactant_dict, product_dict,
                     indices_to_species, forward_rate_constants, human='no'):
    """ This function writes the list of rate expressions for each reaction.
    Parameters
    ----------
    reactant_dict      : dict
                         A dictionary of reactants: key is the reaction number, 
                         values are the list of list with reactant and stoic co-eff
    product_dict      : dict
                         A dictionary of productss: key is the reaction number, 
                         values are the list of list with product and stoic co-eff
    indices_to_species : dict
                         the reverse of speciesindices (keys are the indices
                         and values are the species)
    forward_rate_constants: list
                            A list of forward rate constants for all the reactions
    human              : str, optional
                         indicate whether the output of this function should
                         be formatted for a human to read ('yes'). Default
                         is 'no'
    Returns
    -------
    rates_list_forward : list
                        a list of the rate expressions for all the forward
                        reactions in the model
    rates_list_reverse : list
                        a list of the rate expressions for all the reverse
                        reactions in the model
    """
    kmatrix = forward_rate_constants
    rates_list_forward = build_rate_eqn(kmatrix, reactant_dict, human, forward = 'yes')
    rates_list_reverse = build_rate_eqn(kmatrix, product_dict, human, forward = 'no')
    
    return rates_list_forward, rates_list_reverse
