#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Functions to execute freesasa and parse its output.
"""

from __future__ import division, print_function

import os
import subprocess  # nosec
import sys
import tempfile
from importlib.resources import files

try:
    from Bio.PDB import PDBIO, PDBParser, Select
except ImportError as e:
    print(
        "[!] The binding affinity prediction tools require Biopython",
        file=sys.stderr,
    )
    raise ImportError(e)


import contextlib

from .aa_properties import rel_asa


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr
    https://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python/978264#978264
    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """
    oldstdchannel, dest_file = None, None
    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, "w")
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


def execute_freesasa(structure, selection=None):
    """
    Runs the freesasa executable on a PDB file.

    You can get the executable from:
        https://github.com/mittinatten/freesasa

    The binding affinity models are calibrated with the parameter
    set for vdW radii used in NACCESS:
        http://www.ncbi.nlm.nih.gov/pubmed/994183
    """
    io = PDBIO()

    # try to get freesasa paths from environment if not use the
    #  ones defined in config file
    try:
        freesasa, param_f = [
            os.environ[key] for key in ["FREESASA_BIN", "FREESASA_PAR"]
        ]
    except KeyError as err:
        message = (
            "{} not found. In order to use PRODIGY, set the FREESASA_BIN and "
            "FREESASA_PAR environment variables to point to the freesasa "
            "executable and the naccess.config file respectively."
        )
        raise KeyError(message.format(err))

    if not os.path.isfile(freesasa):
        raise IOError(
            "[!] freesasa binary not found at `{0}`".format(freesasa)
        )
    if not os.path.isfile(param_f):
        raise IOError(
            "[!] Atomic radii file not found at `{0}`".format(param_f)
        )

    # Rewrite PDB using Biopython to have a proper format
    # freesasa is very picky with line width (80 characters or fails!)
    # Select chains if necessary
    class ChainSelector(Select):
        def accept_chain(self, chain):
            if selection and chain.id in selection:
                return 1
            elif not selection:
                return 1
            else:
                return 0

    _pdbf = tempfile.NamedTemporaryFile()
    io.set_structure(structure)
    io.save(_pdbf.name, ChainSelector())

    # Run freesasa
    # Save atomic asa output to another temp file
    _outf = tempfile.NamedTemporaryFile()
    cmd = "{0} -o {1} --format=pdb -c {2} {3}".format(
        freesasa, _outf.name, param_f, _pdbf.name
    )
    p = subprocess.Popen(  # nosec
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = p.communicate()

    if p.returncode:
        print("[!] freesasa did not run successfully", file=sys.stderr)
        print(cmd, file=sys.stderr)
        raise Exception(stderr)

    # Rewind & Parse results file
    # Save
    _outf.seek(0)
    asa, rsa = parse_freesasa_output(_outf)

    _pdbf.close()
    _outf.close()

    return asa, rsa


def parse_freesasa_output(fpath):
    """
    Returns per-residue relative accessibility of side-chain and main-chain
    atoms as calculated by freesasa.
    """

    asa_data, rsa_data = {}, {}

    _rsa = rel_asa
    # _bb = {'CA', 'C', 'N', 'O'}

    p = PDBParser(QUIET=1)
    s = p.get_structure("bogus", fpath.name)
    for res in s.get_residues():
        res_id = (res.parent.id, res.resname, res.id[1])
        _, _, total_asa = 0, 0, 0
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

        rsa_data[res_id] = total_asa / _rsa["total"][res.resname]

    return asa_data, rsa_data


def execute_freesasa_api(structure):
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

    asa_data, rsa_data = {}, {}
    _rsa = rel_asa["total"]

    config_path = os.environ.get(
        "FREESASA_PAR",
        str(files("prodigy_prot").joinpath("naccess.config")),
    )
    classifier = Classifier(config_path)

    # classifier = freesasa.Classifier( os.environ["FREESASA_PAR"])
    # Disable
    with stdchannel_redirected(sys.stderr, os.devnull):
        try:
            struct = structureFromBioPDB(
                structure,
                classifier,
            )
            result = calc(struct)
        except AssertionError as e:
            error_message = "" + os.linesep()
            error_message += "[!] Error when running freesasa:" + os.linesep()
            error_message += f"[!] {e}" + os.linesep()
            error_message += (
                "[!] Make sure the atom names in your PDB file match"
                " the canonical naming and belong "
                "to default residues" + os.linesep()
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
