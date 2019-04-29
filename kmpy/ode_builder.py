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
    
    module_dir = os.getcwd().split('ligpy_utils')[0]
    reactionlist_path = module_dir + '/data/complete_reaction_list.dat'
    rateconstantlist_path = module_dir + '/data/complete_rateconstant_list.dat'
    compositionlist_path = module_dir + '/data/compositionlist.dat'

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
