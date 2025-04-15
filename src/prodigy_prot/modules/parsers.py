"""
Functions to read PDB/mmCIF files
"""

import logging
import sys
import warnings
from pathlib import Path
from typing import Optional, Union

from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.PDB.Structure import Structure

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
log = logging.getLogger("Prodigy")


def get_parser(input_f: Path) -> Union[PDBParser, MMCIFParser]:
    if input_f.suffix == ".cif":
        return MMCIFParser()
    else:
        return PDBParser()


def validate_structure(
    s: Structure, selection: Optional[list[str]] = None, clean: bool = True
) -> Structure:

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
    chains: list[Chain] = list(s.get_chains())
    chain_ids = set([c.id for c in chains])

    if selection:
        sel_chains = []
        # Match selected chain with structure
        for sel in selection:
            for c_str in sel.split(","):
                sel_chains.append(c_str)
                if c_str not in chain_ids:
                    raise ValueError(
                        f"Selected chain not present in provided structure: {c_str}"
                    )

        # Remove unselected chains
        def _ignore_helper(x) -> bool:
            return x.id not in sel_chains

        for c in chains:
            if _ignore_helper(c):
                if c.parent is not None:
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
                    "Unsupported non-standard amino acid found: {0}".format(res.resname)
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
                "{2.parent.id} {2.resname}{2.id[1]}\n".format(i_pp, pp[0], pp[-1])
            )
        logger.warning(message)
        # raise Exception(message)

    return s


def parse_structure(path: str) -> tuple[Structure, int, int]:
    """Return a validated `Structure`, number of chains and number of residues"""

    extension = Path(path).suffix
    supported_extensions = [".pdb", ".cif", ".ent"]
    if extension not in supported_extensions:
        log.error(
            f"[!] Structure format '{extension}' is "
            "not supported. Use '.pdb' or '.cif'."
        )
        sys.exit(1)

    parser = get_parser(Path(path))
    structure_name = Path(path).stem
    structure_path = Path(path)
    try:
        original_structure = parser.get_structure(structure_name, structure_path)
    except Exception as e:
        log.exception(e)
        sys.exit(1)

    assert isinstance(original_structure, Structure)

    structure = validate_structure(original_structure)

    # Get number of chains
    number_of_chains = len(set([c.id for c in structure.get_chains()]))

    # Get number of residues
    number_of_residues = len(list(structure.get_residues()))

    # structure, n_chains, n_res = parse_structure(path=str(struct_path))
    return (structure, number_of_chains, number_of_residues)
