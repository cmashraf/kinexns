import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath + '/../')

import unittest
#import sys
#sys.path.append('../kmpy')

import numpy as np
from ..ode_builder import set_paths, build_species_list
from ..constants import GAS_CONST

paths = set_paths(myPath)


class TestSetPaths(unittest.TestCase):
    """Tests for set_paths()"""

    def test_returns_three_paths(self):
        """Does set_paths() return the expected number of paths?"""
        self.assertEqual(3, len(paths))

    def test_paths_end_correctly(self):
        """Are the correct files in the data directory set?"""
        self.assertEqual('complete_reaction_list.dat', paths[0].split('/')[-1])
        self.assertEqual('complete_rateconstant_list.dat',
                         paths[1].split('/')[-1])
        self.assertEqual('compositionlist.dat', paths[2].split('/')[-1])


class TestGetSpecieslist(unittest.TestCase):
    """Tests for get_specieslist()"""

    def test_correct_num_species(self):
        """Does get_specieslist() return the correct number of species (with
        'correct' being defined as the number of species in the model as
        developed by Hough et al., 2016)"""
        specieslist = build_species_list(paths[0])
        self.assertEqual(6, len(specieslist[2]))

    def test_correct_format(self):
        """Are the entries in specieslist strings longer than 1 character?"""
        specieslist = build_species_list(paths[0])
        for i in np.random.randint(0, len(specieslist[2]), 5):
            self.assertIsInstance(specieslist[2][i], str, msg='%s is not a '
                                                             'string' %
                                                             specieslist[2][i])
            self.assertLess(len(specieslist[2][i]), 2)


    def test_key_value_pairs_are_correct(self):
        """Are the key value pairs in both dictionaries what we expect given
        the original kinetic scheme from Hough et al, 2016?"""
        specieslist = build_species_list(paths[0])
        output = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        self.assertEqual(output['B'], 1)
        self.assertEqual(output['D'], 3)
        self.assertEqual(output['F'], 5)
