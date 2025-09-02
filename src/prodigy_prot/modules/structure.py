from dataclasses import dataclass

from Bio.PDB.Structure import Structure


@dataclass
class ProdigyInput:
    """Class that represents the input to Prodigy."""

    structure: Structure

    def predict(self):
        pass
