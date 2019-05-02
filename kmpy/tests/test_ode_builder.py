import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath + '/../')

import unittest
#import sys
#sys.path.append('../kmpy')

import numpy as np
from ..ode_builder import set_paths, build_species_list
from ..ode_builder import build_kmatrix_forward, build_reac_prod_dict
from ..ode_builder import build_reac_species_dict, build_rate_eqn
from ..ode_builder import build_dydt_list
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

class TestBuildReactantDict(unittest.TestCase):
    """Tests for build_reactant_dict"""

    def test_correct_num_keys(self):
        """Does the returned reactant_dict have the correct number of keys?"""
        specieslist = build_species_list(paths[0])
        output_dict = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        reactantdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        self.assertEqual(len(reactantdict[0]), 4)

    def test_returns_expected_values(self):
        """Does the reactant_dict have the values we expect?"""
        specieslist = build_species_list(paths[0])
        output_dict = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        reactantdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        self.assertEqual(reactantdict[0][1], [[0, 1.0]])
        self.assertEqual(reactantdict[0][2], [[1, 1.0], [2, 1.0]])
        self.assertEqual(reactantdict[1][1], [[4, 1.0]])
        self.assertEqual(reactantdict[1][3], [[3, 1.0], [5, 1.0]])


class TestBuildSpeciesRxnsDict(unittest.TestCase):
    """Tests for build_species_rxns_dict()"""

    def test_correct_num_keys(self):
        """Does the returned reactant_dict have the correct number of keys?"""
        specieslist = build_species_list(paths[0])
        reac_prod_list = [react + prod for react, prod in zip(specieslist[0], specieslist[1])]
        species_rxns = build_reac_species_dict(reac_prod_list, specieslist[2])
        self.assertEqual(len(species_rxns), 6)

    def test_returns_expected_values(self):
        """Does the reactant_dict have the values we expect?"""
        specieslist = build_species_list(paths[0])
        reac_prod_list = [react + prod for react, prod in zip(specieslist[0], specieslist[1])]
        species_rxns = build_reac_species_dict(reac_prod_list, specieslist[2])
        self.assertEqual(species_rxns['B'], [[0, -1, '-1', '+1'], [2, -1, '-1', '+1']])
        self.assertEqual(species_rxns['D'], [[0, 1, '+1', '-1'], [2, 1, '+1', '-1'], [3, 1, '+1', '-1']])
        self.assertEqual(species_rxns['F'], [[3, 1, '+1', '-1']])


class TestBuildRatesEqn(unittest.TestCase):
    """Tests for build_rates_list()"""

    
    def test_correct_length(self):
        """Are the correct number of entries present in the rates_list?"""
        specieslist = build_species_list(paths[0])
        kmatrix = build_kmatrix_forward(paths[1], 298)
        output_dict = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        output_dict_rev = dict(zip(output_dict.values(), output_dict.keys()))
        reacdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        rates_list = build_rate_eqn(kmatrix, reacdict[0], output_dict_rev, human = 'no', forward = 'yes')
        self.assertEqual(len(rates_list), 4)

    def test_returns_expected_values(self):
        """Are the expected rate expressions returned?"""
        kmatrix = build_kmatrix_forward(paths[1], 298)
        specieslist = build_species_list(paths[0])
        output_dict = {specieslist[2][i]:i for i in range(0, len(specieslist[2]))}
        output_dict_rev = dict(zip(output_dict.values(), output_dict.keys()))
        reacdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        rates_list_f = build_rate_eqn(kmatrix, reacdict[0], output_dict_rev, human = 'no', forward = 'yes')
        rates_list_r = build_rate_eqn(kmatrix, reacdict[1], output_dict_rev, human = 'no', forward = 'no')
        self.assertEqual(rates_list_f[0], 'rate_f[0] = kf(T,0) * y[0]**1.0 * y[1]**2.0 ')
        self.assertEqual(rates_list_r[0], 'rate_r[0] = kr(T,0) * y[2]**1.0 * y[3]**1.0 ')


class TestBuildDYDTList(unittest.TestCase):
    """Tests for build_dydt_list()"""

    def test_correct_length(self):
        """Are the correct number of ODEs present in the dydt list?"""
        specieslist = build_species_list(paths[0])
        reac_prod_list = [react + prod for react, prod in zip(specieslist[0], specieslist[1])]
        output_dict = {build_species_list(paths[0])[2][i]:i for i in range(0, len(build_species_list(paths[0])[2]))}
        species_rxns = build_reac_species_dict(reac_prod_list, specieslist[2])
        reacdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        kmatrix = build_kmatrix_forward(paths[1], 298)
        output_dict_rev = dict(zip(output_dict.values(), output_dict.keys()))
        rates_f = build_rate_eqn(kmatrix, reacdict[0], output_dict_rev, human = 'no', forward = 'yes')
        rates_r = build_rate_eqn(kmatrix, reacdict[1], output_dict_rev, human = 'no', forward = 'no')
        dydt = build_dydt_list(rates_f, rates_r, specieslist[2], species_rxns, human='no')    
        self.assertEqual(len(dydt), 6)

    def test_returns_expected_values(self):
        specieslist = build_species_list(paths[0])
        reac_prod_list = [react + prod for react, prod in zip(specieslist[0], specieslist[1])]
        output_dict = {build_species_list(paths[0])[2][i]:i for i in range(0, len(build_species_list(paths[0])[2]))}
        species_rxns = build_reac_species_dict(reac_prod_list, specieslist[2])
        reacdict = build_reac_prod_dict(specieslist[0], specieslist[1], output_dict)
        kmatrix = build_kmatrix_forward(paths[1], 298)
        output_dict_rev = dict(zip(output_dict.values(), output_dict.keys()))
        rates_f = build_rate_eqn(kmatrix, reacdict[0], output_dict_rev, human = 'no', forward = 'yes')
        rates_r = build_rate_eqn(kmatrix, reacdict[1], output_dict_rev, human = 'no', forward = 'no')
        dydt = build_dydt_list(rates_f, rates_r, specieslist[2], species_rxns, human='no')    
        self.assertEqual(dydt[5], 'd[F]/dt = +1*kf(T,3) * y[0]**1.0 -1*kr(T,3) * y[3]**1.0 * y[5]**1.0 ')
