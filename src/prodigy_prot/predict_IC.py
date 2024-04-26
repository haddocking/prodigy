# coding: utf-8
# !/usr/bin/env python
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

from __future__ import division, print_function


__author__ = ["Anna Vangone", "Joao Rodrigues", "Joerg Schaarschmidt"]

import logging
import sys

try:
    from Bio.PDB import NeighborSearch
except ImportError as e:
    print(
        "[!] The binding affinity prediction tools require Biopython",
        file=sys.stderr,
    )
    raise ImportError(e)

from prodigy_prot.modules import aa_properties
from prodigy_prot.modules.freesasa_tools import execute_freesasa_api
from prodigy_prot.modules.models import IC_NIS
from prodigy_prot.modules.parsers import parse_structure
from prodigy_prot.modules.utils import check_path, dg_to_kd


def calculate_ic(struct, d_cutoff=5.5, selection=None):
    """
    Calculates intermolecular contacts in a parsed struct object.
    """
    atom_list = list(struct.get_atoms())
    ns = NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level="R")

    if selection:
        _sd = selection

        def _chain(x):
            return x.parent.id

        ic_list = [
            c
            for c in all_list
            if (_chain(c[0]) in _sd and _chain(c[1]) in _sd)
            and (_sd[_chain(c[0])] != _sd[_chain(c[1])])
        ]
    else:
        ic_list = [c for c in all_list if c[0].parent.id != c[1].parent.id]

    if not ic_list:
        raise ValueError("No contacts found for selection")

    return ic_list


def analyse_contacts(contact_list):
    """
    Enumerates and classifies contacts based on the chemical characteristics
    of the participating amino acids.
    """

    bins = {
        "AA": 0,
        "PP": 0,
        "CC": 0,
        "AP": 0,
        "CP": 0,
        "AC": 0,
    }

    _data = aa_properties.aa_character_ic
    for res_i, res_j in contact_list:
        contact_type = (_data.get(res_i.resname), _data.get(res_j.resname))
        contact_type = "".join(sorted(contact_type))
        bins[contact_type] += 1

    return bins


def analyse_nis(sasa_dict, acc_threshold=0.05):
    """
    Returns the percentages of apolar, polar, and charged
    residues at the interface, according to an accessibility
    criterion.
    """

    _data = aa_properties.aa_character_protorp

    def _char_to_index(x):
        return {"A": 0, "C": 1, "P": 2}.get(x)

    count = [0, 0, 0]

    for res, rsa in sasa_dict.items():
        chain, resn, resi = res
        if rsa >= acc_threshold:
            aa_character = _data[resn]
            aa_index = _char_to_index(aa_character)
            count[aa_index] += 1

    percentages = [100 * x / sum(count) for x in count]
    # print('[+] No. of buried interface residues: {0}'.format(sum(count)))
    return percentages


class Prodigy:
    # init parameters
    def __init__(self, struct_obj, selection=None, temp=25.0):
        self.temp = float(temp)
        if selection is None:
            self.selection = [chain.id for chain in struct_obj.get_chains()]
        else:
            self.selection = selection
        self.structure = struct_obj
        self.ic_network = {}
        self.bins = {}
        self.nis_a = 0
        self.nis_c = 0
        self.ba_val = 0
        self.kd_val = 0

    def predict(self, temp=None, distance_cutoff=5.5, acc_threshold=0.05):
        if temp is not None:
            self.temp = temp
        # Make selection dict from user option or PDB chains
        selection_dict = {}
        for igroup, group in enumerate(self.selection):
            chains = group.split(",")
            for chain in chains:
                if chain in selection_dict:
                    errmsg = (
                        "Selections must be disjoint sets: "
                        f"{chain} is repeated"
                    )
                    raise ValueError(errmsg)
                selection_dict[chain] = igroup

        # Contacts
        self.ic_network = calculate_ic(
            self.structure, d_cutoff=distance_cutoff, selection=selection_dict
        )

        self.bins = analyse_contacts(self.ic_network)

        # SASA
        _, cmplx_sasa = execute_freesasa_api(self.structure)
        self.nis_a, self.nis_c, _ = analyse_nis(
            cmplx_sasa, acc_threshold=acc_threshold
        )

        # Affinity Calculation
        self.ba_val = IC_NIS(
            self.bins["CC"],
            self.bins["AC"],
            self.bins["PP"],
            self.bins["AP"],
            self.nis_a,
            self.nis_c,
        )
        self.kd_val = dg_to_kd(self.ba_val, self.temp)

    def as_dict(self):
        return_dict = {
            "structure": self.structure.id,
            "selection": self.selection,
            "temp": self.temp,
            "ICs": len(self.ic_network),
            "nis_a": self.nis_a,
            "nis_c": self.nis_c,
            "ba_val": self.ba_val,
            "kd_val": self.kd_val,
        }
        return_dict.update(self.bins)
        return return_dict

    def print_prediction(self, outfile="", quiet=False):
        if outfile:
            handle = open(outfile, "w")
        else:
            handle = sys.stdout

        if quiet:
            handle.write(
                "{0}\t{1:8.3f}\n".format(self.structure.id, self.ba_val)
            )
        else:
            handle.write(
                "[+] No. of intermolecular contacts: {0}\n".format(
                    len(self.ic_network)
                )
            )
            handle.write(
                "[+] No. of charged-charged contacts: {0}\n".format(
                    self.bins["CC"]
                )
            )
            handle.write(
                "[+] No. of charged-polar contacts: {0}\n".format(
                    self.bins["CP"]
                )
            )
            handle.write(
                "[+] No. of charged-apolar contacts: {0}\n".format(
                    self.bins["AC"]
                )
            )
            handle.write(
                "[+] No. of polar-polar contacts: {0}\n".format(
                    self.bins["PP"]
                )
            )
            handle.write(
                "[+] No. of apolar-polar contacts: {0}\n".format(
                    self.bins["AP"]
                )
            )
            handle.write(
                "[+] No. of apolar-apolar contacts: {0}\n".format(
                    self.bins["AA"]
                )
            )
            handle.write(
                "[+] Percentage of apolar NIS residues: {0:3.2f}\n".format(
                    self.nis_a
                )
            )
            handle.write(
                "[+] Percentage of charged NIS residues: {0:3.2f}\n".format(
                    self.nis_c
                )
            )
            handle.write(
                "[++] Predicted binding "
                "affinity (kcal.mol-1): {0:8.1f}\n".format(self.ba_val)
            )
            handle.write(
                "[++] Predicted dissociation constant (M) at {:.1f}˚C:"
                " {:8.1e}\n".format(self.temp, self.kd_val)
            )

        if handle is not sys.stdout:
            handle.close()

    def print_contacts(self, outfile=""):
        if outfile:
            handle = open(outfile, "w")
        else:
            handle = sys.stdout

        for res1, res2 in self.ic_network:
            _fmt_str = (
                "{0.resname:>5s} {0.id[1]:5} {0.parent.id:>3s} {1.resname:>5s}"
                " {1.id[1]:5} {1.parent.id:>3s}\n"
            )
            if res1.parent.id not in self.selection[0]:
                res1, res2 = res2, res1
            handle.write(_fmt_str.format(res1, res2))

        if handle is not sys.stdout:
            handle.close()

    def print_pymol_script(self, outfile=""):
        # Writing output PYMOL: pml script
        # initialize array with chains and save chain selection string
        selection_strings = []
        chains = {}
        for s in self.selection:
            selection_strings.append(s.replace(",", "+"))
            for c in s.split(","):
                chains[c] = set()

        # loop over pairs and add interface residues to respective chains
        for pair in self.ic_network:
            for r in pair:
                chains[r.parent.id].add(str(r.id[1]))

        # set output stream
        handle = open(outfile, "w") if outfile else sys.stdout

        # write default setup strings
        handle.writelines(
            [
                "color silver\n",
                "as cartoon\n",
                "bg_color white\n",
                "center\n",
                "color lightblue, chain {}\n".format(selection_strings[0]),
                "color lightpink, chain {}\n".format(selection_strings[1]),
            ]
        )

        # loop over interfaces construct selection strings
        #  and write interface related commands
        for color, iface in [("blue", 1), ("hotpink", 2)]:
            p_sel_string = " or ".join(
                [
                    "chain {} and resi {}".format(c, "+".join(chains[c]))
                    for c in selection_strings[iface - 1].split("+")
                ]
            )
            handle.write("select iface{},  {}\n".format(iface, p_sel_string))
            handle.write("color {}, iface{}\n".format(color, iface))
            handle.write("show sticks, iface{}\n".format(iface))

        # close file handle if applicable
        if handle is not sys.stdout:
            handle.close()


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
    ap.add_argument(
        "structf", help="Structure to analyse in PDB or mmCIF format"
    )
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
    logging.basicConfig(
        level=log_level, stream=sys.stdout, format="%(message)s"
    )
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
