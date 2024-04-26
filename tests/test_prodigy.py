# !/usr/bin/env python

"""
Run the prodigy tests.
"""

from __future__ import division, print_function

import json
import tarfile
import unittest
from io import BufferedReader, TextIOWrapper
from os import devnull
from os.path import basename, dirname, join, splitext
from sys import stderr, version_info

import numpy as np
from Bio.PDB import PDBParser, Structure

from prodigy.modules.freesasa_tools import stdchannel_redirected
from prodigy.modules.parsers import validate_structure
from prodigy.predict_IC import Prodigy


def get_data_path(path):
    """
    Get the path of a file in data for file input

    Args:
      path: path within data folder

    Returns:
        The path to the data file in unix notation
    """
    try:
        return join(dirname(__file__), "..", "data", *path.split("/"))
    except NameError:
        #
        return join("data", *path.split("/"))


class ProdigyOutputTest(unittest.TestCase):
    def test_dataset(self):
        """
        Test method to compare prediction for 80 dataset cases with
            expected values.
        """
        # load expected data from json
        with open(get_data_path("dataset.json")) as fh:
            expected_data = json.load(fh)

        # load dataset PDBs
        dataset = tarfile.open(get_data_path("dataset.tgz"))
        parser = PDBParser(QUIET=True)

        keys_equal = ["AA", "PP", "CC", "AP", "CP", "AC"]
        diffs = {"ba_val": [], "nis_a": [], "nis_c": []}

        # run prodigy for each dataset in the PDB
        for entry in dataset:
            s_name, s_ext = splitext(basename(entry.name))
            # skip system files in archive
            if not s_name.isalnum() or s_ext != ".pdb":
                continue
            # chains = expected_data[s_name]["Interacting_chains"]
            handle = dataset.extractfile(entry)
            # Wrap filehandle to ensure string file handle in Python 3
            if version_info[0] >= 3:
                handle = TextIOWrapper(BufferedReader(handle)) # type: ignore
            # Suppress gap warnings when parsing structure
            with stdchannel_redirected(stderr, devnull):
                s = validate_structure(
                    parser.get_structure(s_name, handle), selection=["A", "B"]
                )
            # Test for structure object
            self.assertIsInstance(s, Structure.Structure)
            # instantiate Prdigy object,
            #  run prediction and retrieve result dict
            prod = Prodigy(s, selection=["A", "B"])
            prod.predict()
            results = prod.as_dict()
            # check for equality of prdicted interface residues
            for k in keys_equal:
                self.assertEqual(results[k], expected_data[s_name][k])
            # check that NIS and binding afinity values are within 2% of
            #  expected values and add diffs for summary
            for k in diffs.keys():
                delta = abs(results[k] / expected_data[s_name][k] - 1)
                # assume a difference of less then 2%
                self.assertAlmostEqual(delta, 0, delta=0.02)
                diffs[k].append(delta)

        # summarize and print calculated delta values
        print("Deltas:")
        for k in diffs.keys():
            print("  {} max: {:.2%}".format(k, max(diffs[k])))
            print("  {} mean: {:.2%}".format(k, np.mean(diffs[k])))


if __name__ == "__main__":
    unittest.main()
