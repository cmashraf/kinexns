#import sys, os
#myPath = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, myPath + '/../')

import unittest

import numpy as np
from kmpy import ode_builder
from kmpy.constants import GAS_CONST

paths = ode_builder.set_paths()


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


