"""
Functions to read PDB/mmCIF files
"""

import logging
import sys
import typing
import warnings
from pathlib import Path
from typing import Optional, Union

from Bio.PDB.Atom import DisorderedAtom
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Model import Model
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


def ignore(r):
    return r.id[0][0] == "W" or r.id[0][0] == "H"


def validate_structure(
    input_strcture_obj: Structure,
    selection: Optional[list[str]] = None,
    clean: bool = True,
) -> list[Model]:

    result: list[Model] = []
    for model in [m for m in input_strcture_obj.child_list]:

        # process selected chains
        chains: list[Chain] = list(model.get_chains())
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
        for atom in list(model.get_atoms()):
            if atom.is_disordered():
                atom = typing.cast(DisorderedAtom, atom)
                residue = atom.parent
                assert residue is not None
                sel_at = atom.selected_child
                assert sel_at is not None
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
            res_list = list(model.get_residues())

            for res in res_list:
                if ignore(res):
                    chain = res.parent
                    assert chain is not None
                    chain.detach_child(res.id)
                elif not is_aa(res, standard=True):
                    raise ValueError(
                        "Unsupported non-standard amino acid found: {0}".format(
                            res.resname
                        )
                    )

            # Remove Hydrogens
            atom_list = list(model.get_atoms())

            def _ignore(x):
                return x.element == "H"

            for atom in atom_list:
                if _ignore(atom):
                    residue = atom.parent
                    assert residue is not None
                    residue.detach_child(atom.name)

        # Detect gaps and compare with no. of chains
        pep_builder = PPBuilder()
        peptides = pep_builder.build_peptides(model)
        n_peptides = len(peptides)

        if n_peptides != len(chain_ids):
            message = "[!] Structure contains gaps:\n"
            for i_pp, pp in enumerate(peptides):
                message += (
                    "\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > "
                    "{2.parent.id} {2.resname}{2.id[1]}\n".format(i_pp, pp[0], pp[-1])
                )
            log.warning(message)

        result.append(model)

    return result


def parse_structure(path: str) -> tuple[list[Model], int, int]:
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

    models: list[Model] = validate_structure(original_structure)

    # Get number of chains
    chain_dict = {}
    res_dict = {}
    for model in models:
        chain_dict.update({c.id: c for c in model.get_chains()})
        res_dict.update({r.id: r for r in model.get_residues()})

    ## Make sure all models have the same chains
    # Get chain sets for all models
    chain_sets = [set(chain.id for chain in model.get_chains()) for model in models]

    # Check if all sets are identical
    if not all(chain_set == chain_sets[0] for chain_set in chain_sets):
        raise ValueError(
            "Not all models have the same chains. Found chain sets: "
            + ", ".join(str(s) for s in chain_sets)
        )

    res_sets = [set(res.id for res in model.get_residues()) for model in models]

    if not all(res_set == res_sets[0] for res_set in res_sets):
        raise ValueError(
            "Not all models have the same residues. Found residue sets: "
            + ", ".join(str(s) for s in res_sets)
        )

    # structure, n_chains, n_res = parse_structure(path=str(struct_path))
    return (models, len(chain_sets[0]), len(res_sets[0]))
