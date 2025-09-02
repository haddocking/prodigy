import sys
from io import TextIOWrapper
from typing import Optional, TextIO, Union

from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from prodigy_prot.modules import aa_properties
from prodigy_prot.modules.freesasa_tools import execute_freesasa_api
from prodigy_prot.modules.models import IC_NIS
from prodigy_prot.modules.utils import dg_to_kd


def calculate_ic(
    model: Model, d_cutoff: float = 5.5, selection: Optional[dict[str, int]] = None
) -> list:
    """
    Calculates intermolecular contacts in a parsed struct object.
    """
    atom_list = list(model.get_atoms())
    ns = NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level="R")

    assert all_list is not None

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

    ic_list.sort()
    return ic_list


def analyse_contacts(contact_list: list) -> dict[str, float]:
    """
    Enumerates and classifies contacts based on the chemical characteristics
    of the participating amino acids.
    """

    bins = {
        "AA": 0.0,
        "PP": 0.0,
        "CC": 0.0,
        "AP": 0.0,
        "CP": 0.0,
        "AC": 0.0,
    }

    _data = aa_properties.aa_character_ic
    for res_i, res_j in contact_list:
        i = _data.get(res_i.resname)
        j = _data.get(res_j.resname)
        if i is not None and j is not None:
            contact_type = "".join(sorted((i, j)))
            bins[contact_type] += 1

    return bins


def analyse_nis(sasa_dict: dict, acc_threshold: float = 0.05) -> list[float]:
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
        _, resn, _ = res
        if rsa >= acc_threshold:
            aa_character = _data[resn]
            aa_index = _char_to_index(aa_character)
            assert aa_index is not None
            count[aa_index] += 1

    percentages = [100.0 * x / sum(count) for x in count]
    # print('[+] No. of buried interface residues: {0}'.format(sum(count)))
    return percentages


class Prodigy:
    # init parameters
    def __init__(
        self,
        model: Model,
        name: str = "",
        selection: Optional[list[str]] = None,
        temp: float = 25.0,
    ):
        self.temp = float(temp)
        if selection is None:
            self.selection = [chain.id for chain in model.get_chains()]
        else:
            self.selection = selection
        self.model = model
        self.name = name
        self.ic_network: list = []
        self.bins: dict[str, float] = {
            "CC": 0.0,
            "CP": 0.0,
            "AC": 0.0,
            "PP": 0.0,
            "AP": 0.0,
            "AA": 0.0,
        }

        self.nis_a = 0.0
        self.nis_c = 0.0
        self.ba_val = 0.0
        self.kd_val = 0.0

    def predict(
        self,
        temp: Optional[float] = None,
        distance_cutoff: float = 5.5,
        acc_threshold: float = 0.05,
    ):
        if temp is not None:
            self.temp = temp
        # Make selection dict from user option or PDB chains
        selection_dict: dict[str, int] = {}
        for igroup, group in enumerate(self.selection):
            chains = group.split(",")
            for chain in chains:
                if chain in selection_dict:
                    errmsg = "Selections must be disjoint sets: " f"{chain} is repeated"
                    raise ValueError(errmsg)
                selection_dict[chain] = igroup

        # Contacts
        self.ic_network = calculate_ic(
            self.model, d_cutoff=distance_cutoff, selection=selection_dict
        )

        self.bins = analyse_contacts(self.ic_network)

        # SASA
        _, cmplx_sasa = execute_freesasa_api(self.model)
        self.nis_a, self.nis_c, _ = analyse_nis(cmplx_sasa, acc_threshold=acc_threshold)

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

    def as_dict(self) -> dict:
        return_dict = {
            "model": self.model.id,
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

    def print_prediction(self, outfile: str = "", quiet: bool = False) -> None:
        handle: Union[TextIOWrapper, TextIO]
        if outfile:
            handle = open(outfile, "w")
        else:
            handle = sys.stdout

        if quiet:
            handle.write("{0}\t{1:8.3f}\n".format(self.name, self.ba_val))
        else:
            handle.write(
                "[+] No. of intermolecular contacts: {0}\n".format(len(self.ic_network))
            )
            handle.write(
                "[+] No. of charged-charged contacts: {0}\n".format(self.bins["CC"])
            )
            handle.write(
                "[+] No. of charged-polar contacts: {0}\n".format(self.bins["CP"])
            )
            handle.write(
                "[+] No. of charged-apolar contacts: {0}\n".format(self.bins["AC"])
            )
            handle.write(
                "[+] No. of polar-polar contacts: {0}\n".format(self.bins["PP"])
            )
            handle.write(
                "[+] No. of apolar-polar contacts: {0}\n".format(self.bins["AP"])
            )
            handle.write(
                "[+] No. of apolar-apolar contacts: {0}\n".format(self.bins["AA"])
            )

            handle.write(
                "[+] Percentage of apolar NIS residues: {0:3.2f}\n".format(self.nis_a)
            )
            handle.write(
                "[+] Percentage of charged NIS residues: {0:3.2f}\n".format(self.nis_c)
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

    def print_contacts(self, outfile: str = "") -> None:
        handle: Union[TextIOWrapper, TextIO]
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

    def print_pymol_script(self, outfile: str = "") -> None:
        # Writing output PYMOL: pml script
        # initialize array with chains and save chain selection string
        selection_strings = []
        chains: dict[str, set] = {}
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
