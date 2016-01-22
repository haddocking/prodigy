#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Binding affinity predictor based on Intermolecular Contacts (ICs).

PL Kastritis , JPGLM Rodrigues, GE Folkers, R Boelens, AMJJ Bonvin
Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface.
Journal of Molecular Biology, 14, 2632-2652 (2014).
"""

from __future__ import print_function, division

__author__ = ["Panagiotis Kastritis", "Joao Rodrigues"]

import os
import sys

from data import aa_properties
from lib.freesasa import execute_freesasa
from lib.models import NIS
from lib.utils import _check_path
from lib.parsers import parse_structure
from lib import aa_properties

def analyse_nis(sasa_dict, acc_threshold=0.05, selection=None):
    """
    Returns the percentages of apolar, polar, and charged
    residues at the interface, according to an accessibility
    criterion.
    """

    _data = aa_properties.aa_character_protorp
    _char_to_index = lambda x: {'A': 0, 'C': 1, 'P': 2}.get(x)
    count = [0, 0, 0]

    for res, rsa in sasa_dict.iteritems():
        chain, resn, resi = res
        if rsa >= acc_threshold:
            aa_character = _data[resn]
            aa_index = _char_to_index(aa_character)
            count[aa_index] += 1

    percentages = map(lambda x: 100*x/sum(count), count)
    # print('[+] No. of buried interface residues: {0}'.format(sum(count)))
    return percentages

def calculate_interface_atoms(cmplx_asa, free_asa, sasa_diff_threshold=1):
    """
    Calculates number of interface atoms in complex based on
    surface area differences between unbound and bound structures.
    """

    n_int_atoms = 0
    for atom, bound_asa in cmplx_asa.iteritems():
        atom_free = free_asa[atom]
        asa_change = atom_free - bound_asa
        if asa_change >= sasa_diff_threshold:
            n_int_atoms += 1
    return n_int_atoms

if __name__ == "__main__":

    try:
        import argparse
        from argparse import RawTextHelpFormatter
    except ImportError as e:
        print('[!] The binding affinity prediction tools require Python 2.7+', file=sys.stderr)
        raise ImportError(e)

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ap.add_argument('structf', help='Structure to analyse in PDB or mmCIF format')
    ap.add_argument('--acc-threshold', type=float, default=0.05, help='Accessibility threshold for BSA analysis')
    ap.add_argument('-q', '--quiet', action='store_true', help='Outputs only the predicted affinity value')

    _co_help = """
    By default, all intermolecular contacts are taken into consideration,
    a molecule being defined as an isolated group of amino acids sharing
    a common chain identifier. In specific cases, for example
    antibody-antigen complexes, some chains should be considered as a
    single molecule.

    Use the --selection option to provide collections of chains that should
    be considered for the calculation. Separate by a space the chains that
    are to be considered _different_ molecules. Use commas to include multiple
    chains as part of a single group:

    --selection A B => Contacts calculated (only) between chains A and B.
    --selection A,B C => Contacts calculated (only) between chains A and C; and B and C.
    --selection A B C => Contacts calculated (only) between chains A and B; B and C; and A and C.
    """
    sel_opt = ap.add_argument_group('Selection Options', description=_co_help)
    sel_opt.add_argument('--selection', nargs='+', metavar=('A B', 'A,B C'))

    cmd = ap.parse_args()

    if cmd.quiet:
        _stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    # Parse structure
    struct_path = _check_path(cmd.structf)
    structure, n_chains, n_res = parse_structure(struct_path)
    print('[+] Parsed PDB file {0} ({1} chains, {2} residues)'.format(structure.id, n_chains, n_res))

    # Make groups from user option or PDB chains
    if cmd.selection:
        group_list = []
        for group in cmd.selection:
            chains = group.split(',')
            for previous in group_list:
                common = chains & previous
                if common:
                    raise ValueError('Selections must be disjoint: {0} is repeated'.format(common))
            group_list.append(chains)
    else:
        group_list = [c.id for c in structure.get_chains()]

    # Complex SASA
    cmplx_asa, cmplx_rsa = execute_freesasa(structure, selection=group_list)
    _, nis_c, nis_p = analyse_nis(cmplx_rsa, acc_threshold=cmd.acc_threshold, selection=group_list)

    # Interface atoms
    free_asa = {}
    for group in group_list:
        group_asa, _ = execute_freesasa(structure, selection=group)
        free_asa.update(group_asa)

    interface_atoms = calculate_interface_atoms(cmplx_asa, free_asa)

    # Affinity Calculation
    ba_val = NIS(nis_c, nis_p, interface_atoms)
    print('[+] Percentage of polar NIS residues: {0:3.2f}'.format(nis_p))
    print('[+] Percentage of charged NIS residues: {0:3.2f}'.format(nis_c))
    print('[+] No. of (buried) interface atoms: {0}'.format(interface_atoms))
    print('[++] Predicted binding affinity: {0:8.3f}'.format(ba_val))

    if cmd.quiet:
        sys.stdout = _stdout
        print('{0}\t{1:8.3f}'.format(struct_path, ba_val))
