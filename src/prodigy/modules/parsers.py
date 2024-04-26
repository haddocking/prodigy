#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Functions to read PDB/mmCIF files
"""

from __future__ import division, print_function

import logging
import os
import sys

try:
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import PPBuilder, is_aa
except ImportError as e:
    print(
        "[!] The binding affinity prediction tools require Biopython",
        file=sys.stderr,
    )
    raise ImportError(e)


def validate_structure(s, selection=None, clean=True):
    # setup logging
    logger = logging.getLogger("Prodigy")

    # Keep first model only
    if len(s) > 1:
        logger.warning(
            (
                "[!] Structure contains more than one model."
                " Only the first one will be kept"
            )
        )
        model_one = s[0].id
        for m in s.child_list[:]:
            if m.id != model_one:
                s.detach_child(m.id)

    # process selected chains
    chains = list(s.get_chains())
    chain_ids = set([c.id for c in chains])

    if selection:
        sel_chains = []
        # Match selected chain with structure
        for sel in selection:
            for c in sel.split(","):
                sel_chains.append(c)
                if c not in chain_ids:
                    raise ValueError(
                        "Selected chain not present"
                        f" in provided structure: {c}"
                    )

        # Remove unselected chains
        def _ignore(x):
            return x.id not in sel_chains

        for c in chains:
            if _ignore(c):
                c.parent.detach_child(c.id)

    # Double occupancy check
    for atom in list(s.get_atoms()):
        if atom.is_disordered():
            residue = atom.parent
            sel_at = atom.selected_child
            sel_at.altloc = " "
            sel_at.disordered_flag = 0
            residue.detach_child(atom.id)
            residue.add(sel_at)

    # Insertion code check
    for c in chains:
        for residue in c.get_residues():
            if residue.get_id()[2] != " ":
                c.detach_child(residue.id)

    if clean:
        # Remove HETATMs and solvent
        res_list = list(s.get_residues())

        def _ignore(r):
            return r.id[0][0] == "W" or r.id[0][0] == "H"

        for res in res_list:
            if _ignore(res):
                chain = res.parent
                chain.detach_child(res.id)
            elif not is_aa(res, standard=True):
                raise ValueError(
                    "Unsupported non-standard amino acid found: {0}".format(
                        res.resname
                    )
                )

        # Remove Hydrogens
        atom_list = list(s.get_atoms())

        def _ignore(x):
            return x.element == "H"

        for atom in atom_list:
            if _ignore(atom):
                residue = atom.parent
                residue.detach_child(atom.name)

    # Detect gaps and compare with no. of chains
    pep_builder = PPBuilder()
    peptides = pep_builder.build_peptides(s)
    n_peptides = len(peptides)

    if n_peptides != len(chain_ids):
        message = "[!] Structure contains gaps:\n"
        for i_pp, pp in enumerate(peptides):
            message += (
                "\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > "
                "{2.parent.id} {2.resname}{2.id[1]}\n".format(
                    i_pp, pp[0], pp[-1]
                )
            )
        logger.warning(message)
        # raise Exception(message)

    return s


def parse_structure(path):
    """
    Parses a structure using Biopython's PDB/mmCIF Parser
    Verifies the integrity of the structure (gaps) and its
    suitability for the calculation (is it a complex?).
    """
    # setup logging
    logger = logging.getLogger("Prodigy")
    logger.info("[+] Reading structure file: {0}".format(path))
    fname = os.path.basename(path)
    sname = ".".join(fname.split(".")[:-1])
    s_ext = fname.split(".")[-1]

    _ext = {"pdb", "ent", "cif"}
    if s_ext not in _ext:
        raise IOError(
            f"[!] Structure format '{s_ext}' is "
            "not supported. Use '.pdb' or '.cif'."
        )

    sparser = PDBParser(QUIET=1) if s_ext in {"pdb", "ent"} else MMCIFParser()

    try:
        s = sparser.get_structure(sname, path)
    except Exception as exeption:
        logger.error("[!] Structure '{0}' could not be parsed".format(sname))
        raise Exception(exeption)

    return (
        validate_structure(s),
        len(set([c.id for c in s.get_chains()])),
        len(list(s.get_residues())),
    )
