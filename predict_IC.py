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

from __future__ import print_function, division

__author__ = ["Anna Vangone", "Joao Rodrigues"]

import os
import subprocess
import sys
import tempfile

try:
    from Bio.PDB import PDBParser, NeighborSearch
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

    # Detect gaps and compare with no. of chains
    pep_builder = PPBuilder()
    peptides = pep_builder.build_peptides(s)
    n_peptides = len(peptides)
    n_chains = len(set([c.id for c in s.get_chains()]))

    if n_peptides != n_chains:
        print('[!] Structure contains gaps:', file=sys.stderr)
        for i_pp, pp in enumerate(peptides):
            print('\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > {2.parent.id} {2.resname}{2.id[1]}'.format(i_pp, pp[0], pp[-1]), file=sys.stderr)
        raise Exception('Calculation cannot proceed')

    return (s, n_chains, n_res)

def calculate_ic(structure, d_cutoff=5.5, selection=None):
    """
    Calculates intermolecular contacts in a parsed structure object.
    """
    atom_list = list(structure.get_atoms())
    ns = NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level='R')

    if selection:
        _sd = selection_dict
        _chain = lambda x: x.parent.id
        ic_list = [c for c in all_list if (_chain(c[0]) in _sd and _chain(c[1]) in _sd)
                    and (_sd[_chain(c[0])] != _sd[_chain(c[1])]) ]
    else:
        ic_list = [c for c in all_list if c[0].parent.id != c[1].parent.id]

    if not ic_list:
        raise ValueError('No contacts found for selection')

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

    _data = aa_properties.aa_character_ic
    for (res_i, res_j) in contact_list:
        contact_type = (_data.get(res_i.resname), _data.get(res_j.resname))
        contact_type = ''.join(sorted(contact_type))
        bins[contact_type] += 1

    return bins

def parse_freesasa_output(fpath):
    """
    Returns per-residue relative accessibility of side-chain and main-chain
    atoms as calculated by freesasa.
    """

    rsa_data = {}

    _rsa = aa_properties.rel_asa
    _bb = set(('CA', 'C', 'N', 'O'))

    s = P.get_structure('bogus', fpath.name)
    for res in s.get_residues():
        res_id = (res.parent.id, res.resname, res.id[1])
        asa_mc, asa_sc, total_asa = 0, 0, 0
        for atom in res:
            asa = atom.bfactor
            # if atom.name in _bb:
            #     asa_mc += asa
            # else:
            #     asa_sc += asa
            total_asa += asa

        rsa_data[res_id] = total_asa / _rsa['total'][res.resname]

    return rsa_data

def get_freesasa():
    """
    Use FREESASA_SRC to get the path to the executable
    and config file.
    """

    src_dir = os.environ.get('FREESASA_SRC')
    if not src_dir:
        print('[!] FREESASA_SRC undefined; Read the installation instructions.', file=sys.stderr)
        raise ValueError('[!] Undefined environment variable: FREESASA_SRC')

    exe_path = os.path.join(src_dir, 'src', 'freesasa')
    _check_path(exe_path)

    config_path = os.path.join(src_dir, 'share', 'naccess.config')
    _check_path(config_path)

    return exe_path, config_path

def execute_freesasa(structure, pdb_selection=None):
    """
    Runs the freesasa executable on a PDB file.

    You can get the executable from:
        https://github.com/mittinatten/freesasa

    The binding affinity models are calibrated with the parameter
    set for vdW radii used in NACCESS:
        http://www.ncbi.nlm.nih.gov/pubmed/994183
    """

    # freesasa = which('freesasa')
    freesasa, param_f= get_freesasa()

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
    cmd = '{0} -B {1} -L -d 0.05 -c {2} {3}'.format(freesasa, _outf.name, param_f, _pdbf.name)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if p.returncode:
        print('[!] freesasa did not run successfully', file=sys.stderr)
        raise Exception(stderr)

    # Rewind & Parse results file
    # Save
    _outf.seek(0)
    rsa = parse_freesasa_output(_outf)

    _pdbf.close()
    _outf.close()

    return rsa

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

def predict_affinity(ic_cc, ic_ca, ic_pp, ic_pa, p_nis_a, p_nis_c):
    """
    Calculates the predicted binding affinity value
    based on the IC-NIS model.
    """

    return 0.09459*ic_cc + 0.10007*ic_ca + -0.19577*ic_pp + 0.22671*ic_pa \
        + -0.18681*p_nis_a + -0.13810*p_nis_c + 15.9433

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

    P = PDBParser(QUIET=1)
    io = PDBIO()

    # Parse structure
    pdb_path = _check_path(cmd.pdbf)
    structure, n_chains, n_res = parse_structure(pdb_path)
    print('[+] Parsed PDB file {0} ({1} chains, {2} residues)'.format(structure.id, n_chains, n_res))

    # Make selection dict from user option or PDB chains
    if cmd.selection:
        selection_dict = {}
        for igroup, group in enumerate(cmd.selection):
            chains = group.split(',')
            for chain in chains:
                if chain in selection_dict:
                    errmsg = 'Selections must be disjoint sets: {0} is repeated'.format(chain)
                    raise ValueError(errmsg)
                selection_dict[chain] = igroup
    else:
        selection_dict = dict([(c.id, nc) for nc, c in enumerate(structure.get_chains())])

    # Contacts
    ic_network = calculate_ic(structure, d_cutoff=cmd.distance_cutoff, selection=selection_dict)
    print('[+] No. of intermolecular contacts: {0}'.format(len(ic_network)))

    bins = analyse_contacts(ic_network)

    # SASA
    cmplx_sasa = execute_freesasa(structure, pdb_selection=selection_dict)
    nis_a, nis_c, _ = analyse_nis(cmplx_sasa, acc_threshold=cmd.acc_threshold, selection=selection_dict)

    # Affinity Calculation
    ba_val = predict_affinity(bins['CC'], bins['AC'], bins['PP'], bins['AP'], nis_a, nis_c)
    print('[+] No. of charged-charged contacts: {0}'.format(bins['CC']))
    print('[+] No. of apolar-charged contacts: {0}'.format(bins['AC']))
    print('[+] No. of polar-polar contacts: {0}'.format(bins['PP']))
    print('[+] No. of apolar-polar contacts: {0}'.format(bins['AP']))
    print('[+] Percentage of apolar NIS residues: {0:3.2f}'.format(nis_a))
    print('[+] Percentage of charged NIS residues: {0:3.2f}'.format(nis_c))
    print('[++] Predicted binding affinity: {0:8.3f}'.format(ba_val))

    if cmd.quiet:
        sys.stdout = _stdout
        print('{0}\t{1:8.3f}'.format(pdb_path, ba_val))
