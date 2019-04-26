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
