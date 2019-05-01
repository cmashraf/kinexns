import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath + '/../')

import unittest
#import sys
#sys.path.append('../kmpy')

import numpy as np
from ..ode_builder import set_paths, build_species_list
from ..ode_builder import build_kmatrix_forward
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

    def test_indicestospecies_keys_match_values_in_speciesindices(self):
        """Do the keys in the dictionary indices_to_species match the
        corresponding values in the dictionary speciesindices?"""
        specieslist = build_species_list(paths[0])
        output_dict = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        output_dict_reverse = dict(zip(output_dict.values(), output_dict.keys()))
        for i in range(len(specieslist[2])):
            spec = output_dict_reverse[i]
            self.assertEqual(output_dict[spec], i)

class TestBuildKMatrix(unittest.TestCase):
    """Tests for build_k_matrix() assuming the model as published by Hough
    et al, 2016"""

    def test_correct_length(self):
        """Does build_k_matrix() return a list of the correct length?"""
        self.assertEqual(len(build_kmatrix_forward(paths[1], 298)), 4)

    def test_correct_k_values(self):
        """Are the entries from build_k_matrix() what we expect?"""
        kmatrix = build_kmatrix_forward(paths[1], 298)
        self.assertEqual(kmatrix[0], 2.287460534686544e-09)
        self.assertEqual(kmatrix[1], 1.526499549219028e-28)
        self.assertEqual(kmatrix[2], 7.640134676963482e-48)
        self.assertEqual(kmatrix[3], 3.399012999643825e-67)
