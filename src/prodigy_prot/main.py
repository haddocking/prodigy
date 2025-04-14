"""
Binding affinity predictor based on Intermolecular Contacts (ICs).
"""

from __future__ import division, print_function

import logging
import sys

from prodigy_prot.parsers import parse_structure
from prodigy_prot.utils import check_path
from prodigy_prot.prodigy import Prodigy


def main():
    try:
        import argparse
        from argparse import RawTextHelpFormatter
    except ImportError as err:
        print(
            "[!] The binding affinity prediction tools require Python 2.7+",
            file=sys.stderr,
        )
        raise ImportError(err)

    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    ap.add_argument("structf", help="Structure to analyse in PDB or mmCIF format")
    ap.add_argument(
        "--distance-cutoff",
        type=float,
        default=5.5,
        help="Distance cutoff to calculate ICs",
    )
    ap.add_argument(
        "--acc-threshold",
        type=float,
        default=0.05,
        help="Accessibility threshold for BSA analysis",
    )
    ap.add_argument(
        "--temperature",
        type=float,
        default=25.0,
        help="Temperature (C) for Kd prediction",
    )
    ap.add_argument(
        "--contact_list", action="store_true", help="Output a list of contacts"
    )
    ap.add_argument(
        "--pymol_selection",
        action="store_true",
        help="Output a script to highlight the interface (pymol)",
    )
    ap.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Outputs only the predicted affinity value",
    )
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
    --selection A,B C => Contacts calculated (only) between \
        chains A and C; and B and C.
    --selection A B C => Contacts calculated (only) between \
        chains A and B; B and C; and A and C.
    """
    sel_opt = ap.add_argument_group("Selection Options", description=_co_help)
    sel_opt.add_argument("--selection", nargs="+", metavar=("A B", "A,B C"))

    cmd = ap.parse_args()

    # setup logging
    log_level = logging.ERROR if cmd.quiet else logging.INFO
    logging.basicConfig(level=log_level, stream=sys.stdout, format="%(message)s")
    logger = logging.getLogger("Prodigy")

    struct_path = check_path(cmd.structf)

    # Parse structure
    structure, n_chains, n_res = parse_structure(struct_path)
    logger.info(
        "[+] Parsed structure file {0} ({1} chains, {2} residues)".format(
            structure.id, n_chains, n_res
        )
    )
    prodigy = Prodigy(structure, cmd.selection, cmd.temperature)
    prodigy.predict(
        distance_cutoff=cmd.distance_cutoff, acc_threshold=cmd.acc_threshold
    )
    prodigy.print_prediction(quiet=cmd.quiet)

    # Print out interaction network
    if cmd.contact_list:
        fname = struct_path[:-4] + ".ic"
        prodigy.print_contacts(fname)

    # Print out interaction network
    if cmd.pymol_selection:
        fname = struct_path[:-4] + ".pml"
        prodigy.print_pymol_script(fname)


if __name__ == "__main__":
    main()
