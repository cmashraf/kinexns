import os
# sys.path.insert(0, myPath + '/../')

import unittest
# import sys
# sys.path.append('../kmpy')

import numpy as np
from ..ode_builder import set_paths, build_species_list
from ..ode_builder import build_kmatrix_forward
from ..ode_builder import build_dydt_list
from ..ode_builder import build_free_energy_dict, build_free_energy_change
from ..ode_builder import build_kmatrix_reverse

myPath = os.path.dirname(os.path.abspath(__file__))
# from ..constants import GAS_CONST

paths = set_paths(myPath)
specieslist = build_species_list(paths[0])
reac_prod_list = [react + prod for react, prod in
                  zip(specieslist[0], specieslist[1])]
output_dict = {build_species_list(paths[0])[2][i]: i
               for i in range(0, len(build_species_list(paths[0])[2]))}
kmatrix_f = build_kmatrix_forward(paths[1], 573)
output_dict_rev = dict(zip(output_dict.values(), output_dict.keys()))
free_energy_dict = build_free_energy_dict(paths[2], 573)
free_energy_change, mol_change = build_free_energy_change(reac_prod_list,
                                                          free_energy_dict)
kmatrix_r = build_kmatrix_reverse(reac_prod_list, free_energy_dict,
                                  kmatrix_f, 573)


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
        self.assertEqual('free_energy_library.dat', paths[2].split('/')[-1])


class TestBuildSpecieslist(unittest.TestCase):
    """Tests for build_species_list()"""

    def test_correct_num_species(self):
        """Does get_specieslist() return the correct number of species (with
        'correct' being defined as the number of species in the model as
        developed by Hough et al., 2016)"""
        self.assertEqual(6, len(specieslist[2]))

    def test_correct_format(self):
        """Are the entries in specieslist strings longer than 1 character?"""
        for i in np.random.randint(0, len(specieslist[2]), 5):
            self.assertIsInstance(specieslist[2][i], str, msg=' %s is not'
                                  'a string' % specieslist[2][i])
            self.assertLess(len(specieslist[2][i]), 2)

    def test_key_value_pairs_are_correct(self):
        """Are the key value pairs in both dictionaries what we expect given
        the original kinetic scheme from Hough et al, 2016?"""
        self.assertEqual(output_dict['B'], 1)
        self.assertEqual(output_dict['D'], 3)
        self.assertEqual(output_dict['F'], 5)

    def test_indicestospecies_keys_match_values_in_speciesindices(self):
        """Do the keys in the dictionary indices_to_species match the
        corresponding values in the dictionary speciesindices?"""
        for i in range(len(specieslist[2])):
            spec = output_dict_rev[i]
            self.assertEqual(output_dict[spec], i)


class TestBuildKMatrix(unittest.TestCase):
    """Tests for build_k_matrix() assuming the model as published by Hough
    et al, 2016"""

    def test_correct_length_f(self):
        """Does build_k_matrix() return a list of the correct length?"""
        self.assertEqual(len(build_kmatrix_forward(paths[1], 573)), 4)

    def test_correct_k_values_f(self):
        """Are the entries from build_k_matrix() what we expect?"""
        # kmatrix=build_kmatrix_forward(paths[1], 298)
        self.assertEqual(kmatrix_f[0], 9.070586400466592e-07)
        self.assertEqual(kmatrix_f[1], 9.070586400466592e-07)
        self.assertEqual(kmatrix_f[2], 1.8791556907858517e-07)
        self.assertEqual(kmatrix_f[3], 3.548590208923674e-05)


class TestBuildKMatrixRev(unittest.TestCase):
    """Tests for build_k_matrix() assuming the model as published by Hough
    et al, 2016"""

    def test_correct_length(self):
        """Does build_k_matrix() return a list of the correct length?"""
        self.assertEqual(len(free_energy_dict), len(specieslist[2]))

    def test_correct_length_kmatrix_r(self):
        """Does build_k_matriix() return a list of the correct length?"""
        self.assertEqual(len(kmatrix_r), 4)

    def test_correct_length_free_energy(self):
        """Does build_k_matriix() return a list of the correct length?"""
        self.assertEqual(len(free_energy_change), 4)

    def test_correct_free_energy_values(self):
        self.assertEqual(free_energy_dict['B'], -76.44006077)
        self.assertEqual(free_energy_dict['D'], -496.14661846999996)
        self.assertEqual(free_energy_dict['F'], -572.5674811700001)

    def test_correct_free_energy_cahnge(self):
        self.assertEqual(free_energy_change[0], -6.293953619833474)
        self.assertEqual(free_energy_change[2], -32.32967185957523)

    def test_correct_k_values_r(self):
        """Are the entries from build_k_matrix() what we expect?"""
        # kmatrix=build_kmatrix_forward(paths[1], 298)
        self.assertEqual(kmatrix_r[0], 1.1379057686786213e-05)
        self.assertEqual(kmatrix_r[1], 2.895704050369846e-08)
        self.assertEqual(kmatrix_r[2], 9.975478106112853e-09)
        self.assertEqual(kmatrix_r[3], 0.9483448135351068)
