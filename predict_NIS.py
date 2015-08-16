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
import subprocess
import sys
import tempfile

try:
    from Bio.PDB import PDBParser
    from Bio.PDB import PDBIO, Select
    from Bio.PDB.Polypeptide import PPBuilder, is_aa
except ImportError as e:
    print('[!] The binding affinity prediction tools require Biopython', file=sys.stderr)
    raise ImportError(e)

from data import aa_properties

#
def parse_structure(path):
    """
    Parses a PDB formatter structure using Biopython's PDB Parser
    Verifies the integrity of the structure (gaps) and its
    suitability for the calculation (is it a complex?).
    """

    print('[+] Reading structure file: {0}'.format(path))
    fname = os.path.basename(path)
    sname = '.'.join(fname.split('.')[:-1])

    try:
        s = P.get_structure(sname, path)
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

    nn_res = len(list(s.get_residues()))

    print('[+] Parsed PDB file {0} ({1}/{2} residues kept)'.format(sname, nn_res, n_res))

    # Detect gaps and compare with no. of chains
    pep_builder = PPBuilder()
    peptides = pep_builder.build_peptides(s)
    n_peptides = len(peptides)
    n_chains = len(set([c.id for c in s.get_chains()]))

    if n_peptides != n_chains:
        print('[!] Structure contains gaps:', file=sys.stderr)
        for i_pp, pp in enumerate(peptides):
            print('\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > {2.parent.id} {2.resname}{2.id[1]}'.format(i_pp, pp[0], pp[-1]))
        raise Exception('Calculation cannot proceed')

    return s

def parse_freesasa_output(fpath):
    """
    Returns per-residue relative accessibility of side-chain and main-chain
    atoms as calculated by freesasa.
    """

    asa_data, rsa_data = {}, {}

    _rsa = aa_properties.rel_asa
    _bb = set(('CA', 'C', 'N', 'O'))

    s = P.get_structure('bogus', fpath.name)
    for res in s.get_residues():
        res_id = (res.parent.id, res.resname, res.id[1])
        asa_mc, asa_sc, total_asa = 0, 0, 0
        for atom in res:
            aname = atom.name
            at_id = (res.parent.id, res.resname, res.id[1], aname)
            asa = atom.bfactor
            # if atom.name in _bb:
            #     asa_mc += asa
            # else:
            #     asa_sc += asa
            total_asa += asa
            asa_data[at_id] = asa

        rsa_data[res_id] = total_asa / _rsa['total'][res.resname]

    return asa_data, rsa_data

def execute_freesasa(structure, pdb_selection=None):
    """
    Runs the freesasa executable on a PDB file.

    You can get the executable from:
        https://github.com/mittinatten/freesasa

    The binding affinity models are calibrated with the parameter
    set for vdW radii used in NACCESS:
        http://www.ncbi.nlm.nih.gov/pubmed/994183
    """

    try:
        with open(os.devnull, 'w') as void:
            subprocess.call('freesasa -h', shell=True, stderr=void)
    except OSError as e:
        raise Exception('Could not execute `freesasa`: {0}'.format(e))

    # Rewrite PDB using Biopython to have a proper format
    # freesasa is very picky with line width (80 characters or fails!)
    # Select chains if necessary
    class ChainSelector(Select):
        def accept_chain(self, chain):
            if pdb_selection and chain.id in pdb_selection:
                return 1
            elif not pdb_selection:
                return 1
            else:
                return 0

    _pdbf = tempfile.NamedTemporaryFile()
    io.set_structure(structure)
    io.save(_pdbf.name, ChainSelector())

    # Run freesasa
    # Save atomic asa output to another temp file
    _outf = tempfile.NamedTemporaryFile()
    cmd = 'freesasa -B {0} -L -d 0.05 {1}'.format(_outf.name, _pdbf.name)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    # Rewind & Parse results file
    # Save
    _outf.seek(0)
    asa, rsa = parse_freesasa_output(_outf)

    _pdbf.close()
    _outf.close()

    return asa, rsa

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
    print('[+] No. of buried interface residues: {0}'.format(sum(count)))
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

def predict_affinity(p_nis_c, p_nis_p, n_int_atoms):
    """
    Calculates the predicted binding affinity value
    based on the NIS model.
    """

    return 0.0856851248873*p_nis_p + -0.0685254498746*p_nis_c + 0.0261591389985*n_int_atoms \
        + 3.0124939659498

def _check_path(path):
    """
    Checks if a file is readable.
    """

    full_path = os.path.abspath(path)
    if not os.path.isfile(full_path):
        raise IOError('Could not read file: {0}'.format(path))
    return full_path

if __name__ == "__main__":

    try:
        import argparse
        from argparse import RawTextHelpFormatter
    except ImportError as e:
        print('[!] The binding affinity prediction tools require Python 2.7+', file=sys.stderr)
        raise ImportError(e)

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ap.add_argument('pdbf', help='Structure to analyse in PDB format')
    ap.add_argument('--distance-cutoff', default=5.5, help='Distance cutoff to calculate ICs')
    ap.add_argument('--acc-threshold', default=0.05, help='Accessibility threshold for BSA analysis')
    ap.add_argument('--quiet', action='store_true', help='Outputs only the predicted affinity value')
    # ap.add_argument('--outfile', default=sys.stdout, help='Output file where to write analysis')
    # ap.add_argument('--outfmt', default='tabular', choices=['csv', 'tabular'],
    #                 help='Output file format')

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

    P = PDBParser(QUIET=1)
    io = PDBIO()

    # Parse structure
    pdb_path = _check_path(cmd.pdbf)
    structure = parse_structure(pdb_path)

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
    cmplx_asa, cmplx_rsa = execute_freesasa(structure, pdb_selection=group_list)
    _, nis_c, nis_p = analyse_nis(cmplx_rsa, acc_threshold=cmd.acc_threshold, selection=group_list)

    # Interface atoms
    free_asa = {}
    for group in group_list:
        group_asa, _ = execute_freesasa(structure, pdb_selection=group)
        free_asa.update(group_asa)

    interface_atoms = calculate_interface_atoms(cmplx_asa, free_asa)

    # Affinity Calculation
    ba_val = predict_affinity(nis_c, nis_p, interface_atoms)
    print('[+] Predicted binding affinity: {0:8.3f}'.format(ba_val))

    if cmd.quiet:
        sys.stdout = _stdout
        print('{0}\t{1:8.3f}'.format(pdb_path, ba_val))
