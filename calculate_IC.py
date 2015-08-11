#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Binding affinity predictor based on Intermolecular Contacts (ICs).

Anna Vangone and Alexandre M.J.J. Bonvin,
Contacts-based prediction of binding affinity in protein-protein complexes.
eLife (2015)
"""

from __future__ import print_function

__author__ = ["Anna Vangone", "Joao Rodrigues"]

import os
import sys

try:
    from Bio.PDB import PDBParser, NeighborSearch
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

    full_path = os.path.abspath(path)
    if not os.path.isfile(full_path):
        raise IOError('[!] Could not read input structure: {0}'.format(path), file=sys.stderr)
    else:
        print('[+] Reading structure file: {0}'.format(path))
        fname = os.path.basename(full_path)
        sname = '.'.join(fname.split('.')[:-1])

    try:
        s = P.get_structure(sname, full_path)
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
            raise ValueError('[!] Unsupported non-standard amino acid found: {0}'.format(res.resname))

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

def calculate_ic(structure, d_cutoff):
    """
    Calculates intermolecular contacts in a parsed structure object.
    """
    atom_list = list(structure.get_atoms())
    ns = NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level='R')
    ic_list = [c for c in all_list if c[0].parent.id != c[1].parent.id]

    print('[+] No. of intermolecular contacts: {0}/{1}'.format(len(ic_list), len(all_list)))
    return ic_list

def analyse_contacts(contact_list):
    """
    Enumerates and classifies contacts based on the chemical characteristics
    of the participating amino acids.
    """

    bins = {
        'AA': 0, 'PP': 0,
        'CC': 0, 'AP': 0,
        'CP': 0, 'AC': 0,
        }

    _data = aa_properties.aa_character
    for (res_i, res_j) in contact_list:
        contact_type = (_data.get(res_i.resname), _data.get(res_j.resname))
        contact_type = ''.join(sorted(contact_type))
        bins[contact_type] += 1

    return bins

def write_contacts(contact_data, outfile, oformat):
    """
    Outputs the result of the analysis to a file or stdout.
    """

    # Define format here
    oformat_dict = {
        'csv': '{0},{1}\n',
        'tabular': '{0}\t{1}\n',
        }

    _ostr = oformat_dict[oformat]

    print('[+] Writing analysis output to: {0}'.format(outfile.name))
    for contact_type in sorted(contact_data):
        n_matches = contact_data[contact_type]
        outfile.write(_ostr.format(contact_type, n_matches))

if __name__ == "__main__":

    try:
        import argparse
        from argparse import RawTextHelpFormatter
    except ImportError as e:
        print('[!] The binding affinity prediction tools require Python 2.7+', file=sys.stderr)
        raise ImportError(e)

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ap.add_argument('pdb_list', nargs='+', help='Structure(s) to analyse in PDB format')
    ap.add_argument('--cutoff', default=5.5, help='Distance cutoff to calculate ICs')
    ap.add_argument('--outfile', default=sys.stdout, help='Output file where to write analysis')
    ap.add_argument('--outfmt', default='tabular', choices=['csv', 'tabular'],
                    help='Output file format')

    cmd = ap.parse_args()

    P = PDBParser(QUIET=1)

    for pdbf in cmd.pdb_list:
        structure = parse_structure(pdbf)
        ic_network = calculate_ic(structure, d_cutoff=cmd.cutoff)
        bins = analyse_contacts(ic_network)

        if isinstance(cmd.outfile, str):
            cmd.outfile = open(cmd.outfile, 'w')
        with cmd.outfile as handle:
            write_contacts(bins, handle, cmd.outfmt)
