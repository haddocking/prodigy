# !/usr/bin/env python

"""
Run the prodigy tests.
"""



import json
import unittest
from os.path import join, basename, splitext, dirname
from os import environ,devnull
import tarfile
from io import TextIOWrapper, BufferedReader
from sys import stderr

from Bio.PDB import PDBParser,Structure
from prodigy import Prodigy
from prodigy.predict_IC import main
from prodigy.lib.parsers import validate_structure
from prodigy.lib.freesasa_tools import execute_freesasa, execute_freesasa_api, stdchannel_redirected
import numpy as np

def get_data_path(path):
    try:
        return join(dirname(__file__), "..", "data", *path.split('/'))
    except NameError:
        return join( "data", *path.split('/'))

class ProdigyOutputTest(unittest.TestCase):
    def test_dataset(self):
        # load expected data from json
        with open(get_data_path('dataset.json')) as fh :
            expected_data = json.load(fh)

        # load dataset PDBs
        dataset = tarfile.open(get_data_path('dataset.tgz'))
        parser = PDBParser(QUIET=True)

        keys_equal = ['AA', 'PP', 'CC', 'AP', 'CP', 'AC']
        diffs = {'ba_val':[], 'nis_a':[], 'nis_c':[]}

        # run prodigy for each dataset in the PDB
        for entry in dataset:
            s_name, s_ext = splitext(basename(entry.name))
            # skip system files in archive
            if not s_name.isalnum() or s_ext != '.pdb':
                continue
            # chains = expected_data[s_name]["Interacting_chains"]
            handle = TextIOWrapper(BufferedReader(dataset.extractfile(entry)))
            with stdchannel_redirected(stderr, devnull):
                s = validate_structure(parser.get_structure(s_name, handle),selection=['A','B'])

            self.assertIsInstance(s,Structure.Structure)
            prod = Prodigy(s,selection=['A','B'])
            prod.predict()
            results = prod.as_dict()
            for k in keys_equal:
                self.assertEqual(results[k], expected_data[s_name][k])
            for k in diffs.keys():
                delta = abs(results[k]/expected_data[s_name][k] - 1)
                # assume a difference of less then 2%
                self.assertAlmostEqual(delta, 0, delta=0.015)
                diffs[k].append(delta)
        print('Deltas:')
        for k in diffs.keys():
            print("  {} max: {:.2%}".format(k, max(diffs[k])))
            print("  {} mean: {:.2%}".format(k, np.mean(diffs[k])))


