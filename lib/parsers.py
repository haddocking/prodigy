#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Functions to read PDB/mmCIF files
"""

from __future__ import print_function, division

import os
import sys

try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.Polypeptide import PPBuilder, is_aa
except ImportError as e:
    print('[!] The binding affinity prediction tools require Biopython', file=sys.stderr)
    raise ImportError(e)

def parse_structure(path):
    """
    Parses a structure using Biopython's PDB/mmCIF Parser
    Verifies the integrity of the structure (gaps) and its
    suitability for the calculation (is it a complex?).
    """

    print('[+] Reading structure file: {0}'.format(path))
    fname = os.path.basename(path)
    sname = '.'.join(fname.split('.')[:-1])
    s_ext = fname.split('.')[-1]

    _ext = set(('pdb', 'ent', 'cif'))
    if s_ext not in _ext:
        raise IOError('[!] Structure format \'{0}\' is not supported. Use \'.pdb\' or \'.cif\'.'.format(s_ext))

    if s_ext in set(('pdb', 'ent')):
        sparser = PDBParser(QUIET=1)
    elif s_ext == 'cif':
        sparser = MMCIFParser()

    try:
        s = sparser.get_structure(sname, path)
    except Exception as e:
        print('[!] Structure \'{0}\' could not be parsed'.format(sname), file=sys.stderr)
        raise Exception(e)

    # Double occupancy check
    for atom in list(s.get_atoms()):
        if atom.is_disordered():
            residue = atom.parent
            sel_at = atom.selected_child
            sel_at.altloc = ' '
            sel_at.disordered_flag = 0
            residue.detach_child(atom.id)
            residue.add(sel_at)

    # Remove HETATMs and solvent
    res_list = list(s.get_residues())
    n_res = len(res_list)
    _ignore = lambda r: r.id[0][0] == 'W' or r.id[0][0] == 'H'
    for res in res_list:
        if _ignore(res):
            chain = res.parent
            chain.detach_child(res.id)
        elif not is_aa(res, standard=True):
            raise ValueError('Unsupported non-standard amino acid found: {0}'.format(res.resname))

    # Detect gaps and compare with no. of chains
    pep_builder = PPBuilder()
    peptides = pep_builder.build_peptides(s)
    n_peptides = len(peptides)
    n_chains = len(set([c.id for c in s.get_chains()]))

    if n_peptides != n_chains:
        print('[!] Structure contains gaps:', file=sys.stderr)
        for i_pp, pp in enumerate(peptides):
            print('\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > {2.parent.id} {2.resname}{2.id[1]}'.format(i_pp, pp[0], pp[-1]), file=sys.stderr)
        #raise Exception('Calculation cannot proceed')

    return (s, n_chains, n_res)