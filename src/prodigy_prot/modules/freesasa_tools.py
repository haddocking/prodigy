"""
Functions to execute freesasa and parse its output.
"""

import os

import freesasa
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from freesasa import Classifier, calc, structureFromBioPDB

from prodigy_prot import NACCESS_CONFIG
from prodigy_prot.modules.aa_properties import rel_asa

freesasa.setVerbosity(freesasa.nowarnings)


def execute_freesasa_api(model: Model) -> tuple[dict, dict]:
    """
    Calls freesasa using its Python API and returns
    per-residue accessibilities.
    """

    asa_data = {}
    rsa_data: dict[tuple[str, int, str], float] = {}
    _rsa: dict = rel_asa["total"]

    classifier = Classifier(str(NACCESS_CONFIG))

    # NOTE: `structureFromBioPDB` requires a Structure object
    #  so here build one from a model
    s = Structure(model.id)
    s.add(model)

    try:
        struct = structureFromBioPDB(
            s,
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
