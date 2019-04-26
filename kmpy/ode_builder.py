import os
#import numpy as np
#from constants import GAS_CONST, MW


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
        
        #getting the reactants for each reaction
        for spec in line.split(','):
            if float(spec.split('_')[0].split()[0]) < 0:
                self.reactants_names.append((spec.split('_')[0].split()[0],
                                          spec.split('_')[1].split()[0]))
            #print(self.species_names)
        return self.reactants_names
    
    def getProductsName(self, line):
        
        #getting the reactants for each reaction
        for spec in line.split(','):
            if float(spec.split('_')[0].split()[0]) > 0:
                self.products_names.append((spec.split('_')[0].split()[0],
                                          spec.split('_')[1].split()[0]))
            #print(self.species_names)
        return self.products_names
    
    def uniqueSpeciesName(self, line, species_list):
        
        #building the unique species list
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
