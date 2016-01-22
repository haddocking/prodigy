#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Generic properties of amino acids required for the binding affinity
prediction methods.
"""

__author__ = ["Anna Vangone", "Joao Rodrigues"]

aa_character_ic = {
    'ALA': 'A',
    'CYS': 'A', # ?
    'GLU': 'C',
    'ASP': 'C',
    'GLY': 'A',
    'PHE': 'A',
    'ILE': 'A',
    'HIS': 'C',
    'LYS': 'C',
    'MET': 'A',
    'LEU': 'A',
    'ASN': 'P',
    'GLN': 'P',
    'PRO': 'A',
    'SER': 'P',
    'ARG': 'C',
    'THR': 'P',
    'TRP': 'A',
    'VAL': 'A',
    'TYR': 'A',
}

aa_character_protorp = {
    'ALA': 'A',
    'CYS': 'P',
    'GLU': 'C',
    'ASP': 'C',
    'GLY': 'A',
    'PHE': 'A',
    'ILE': 'A',
    'HIS': 'P',
    'LYS': 'C',
    'MET': 'A',
    'LEU': 'A',
    'ASN': 'P',
    'GLN': 'P',
    'PRO': 'A',
    'SER': 'P',
    'ARG': 'C',
    'THR': 'P',
    'TRP': 'P',
    'VAL': 'A',
    'TYR': 'P',
}

# Scaling factors for relative ASA
# Calculated using extended ALA-X-ALA peptides
# Taken from NACCESS
rel_asa = {
    'total':
        {
        'ALA': 107.95,
        'CYS': 134.28,
        'ASP': 140.39,
        'GLU': 172.25,
        'PHE': 199.48,
        'GLY': 80.10,
        'HIS': 182.88,
        'ILE': 175.12,
        'LYS': 200.81,
        'LEU': 178.63,
        'MET': 194.15,
        'ASN': 143.94,
        'PRO': 136.13,
        'GLN': 178.50,
        'ARG': 238.76,
        'SER': 116.50,
        'THR': 139.27,
        'VAL': 151.44,
        'TRP': 249.36,
        'TYR': 212.76,
        },
    'bb':
        {
        'ALA': 38.54,
        'CYS': 37.53,
        'ASP': 37.70,
        'GLU': 37.51,
        'PHE': 35.37,
        'GLY': 47.77,
        'HIS': 35.80,
        'ILE': 37.16,
        'LYS': 37.51,
        'LEU': 37.51,
        'MET': 37.51,
        'ASN': 37.70,
        'PRO': 16.23,
        'GLN': 37.51,
        'ARG': 37.51,
        'SER': 38.40,
        'THR': 37.57,
        'VAL': 37.16,
        'TRP': 38.10,
        'TYR': 35.38,
        },
    'sc':
        {
        'ALA': 69.41,
        'CYS': 96.75,
        'ASP': 102.69,
        'GLU': 134.74,
        'PHE': 164.11,
        'GLY': 32.33,
        'HIS': 147.08,
        'ILE': 137.96,
        'LYS': 163.30,
        'LEU': 141.12,
        'MET': 156.64,
        'ASN': 106.24,
        'PRO': 119.90,
        'GLN': 140.99,
        'ARG': 201.25,
        'SER': 78.11,
        'THR': 101.70,
        'VAL': 114.28,
        'TRP': 211.26,
        'TYR': 177.38,
    }
}
