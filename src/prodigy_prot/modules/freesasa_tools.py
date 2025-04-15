"""
Functions to execute freesasa and parse its output.
"""

import os
import sys

from Bio.PDB.Structure import Structure

from prodigy_prot import NACCESS_CONFIG
from prodigy_prot.modules.aa_properties import rel_asa


def execute_freesasa_api(structure: Structure) -> tuple[dict, dict]:
    """
    Calls freesasa using its Python API and returns
    per-residue accessibilities.
    """
    try:
        from freesasa import Classifier, calc, structureFromBioPDB
    except ImportError as err:
        print(
            (
                "[!] The binding affinity prediction tools require the "
                "'freesasa' Python API"
            ),
            file=sys.stderr,
        )
        raise ImportError(err)

    asa_data = {}
    rsa_data: dict[tuple[str, int, str], float] = {}
    _rsa: dict = rel_asa["total"]

    classifier = Classifier(str(NACCESS_CONFIG))

    try:
        struct = structureFromBioPDB(
            structure,
            classifier,
        )
        result = calc(struct)
    except AssertionError as e:
        error_message = "" + os.linesep
        error_message += "[!] Error when running freesasa:" + os.linesep
        error_message += f"[!] {e}" + os.linesep
        error_message += (
            "[!] Make sure the atom names in your PDB file match"
            " the canonical naming and belong "
            "to default residues" + os.linesep
        )
        print(error_message)
        raise Exception(error_message)

    # iterate over all atoms to get SASA and residue name
    for idx in range(struct.nAtoms()):
        atname = struct.atomName(idx)
        resname = struct.residueName(idx)
        resid = struct.residueNumber(idx)
        chain = struct.chainLabel(idx)
        at_uid = (chain, resname, resid, atname)
        res_uid = (chain, resname, resid)

        asa = result.atomArea(idx)
        asa_data[at_uid] = asa
        # add asa to residue
        rsa_data[res_uid] = rsa_data.get(res_uid, 0) + asa

    # convert total asa ro relative asa
    rsa_data.update(
        (res_uid, asa / _rsa[res_uid[1]]) for res_uid, asa in rsa_data.items()
    )
    return asa_data, rsa_data
