from pathlib import Path

import pytest
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from prodigy_prot.modules.parsers import parse_structure, validate_structure

from . import TEST_DATA


@pytest.fixture
def input_structure_cif():
    yield Path(TEST_DATA, "2oob.cif")


@pytest.fixture
def input_structure_pdb() -> Path:
    return Path(TEST_DATA, "2oob.pdb")


def test_validate_stucture_pdb(input_structure_pdb):

    parser = PDBParser()
    structure = parser.get_structure("test_structure", input_structure_pdb)
    assert isinstance(structure, Structure)

    result = validate_structure(structure)
    assert result == structure


def test_validate_stucture_cif(input_structure_cif):

    parser = MMCIFParser()
    structure = parser.get_structure("test_structure", input_structure_cif)
    assert isinstance(structure, Structure)

    result = validate_structure(structure)
    assert result == structure


def test_parse_structure_pdb(input_structure_pdb):

    parser = PDBParser()
    structure = parser.get_structure(input_structure_pdb.stem, input_structure_pdb)

    result, num_chains, num_res = parse_structure(input_structure_pdb)

    assert result == structure
    assert num_chains == 2
    assert num_res == 116


def test_parse_structure_cif(input_structure_cif):

    parser = MMCIFParser()
    structure = parser.get_structure(input_structure_cif.stem, input_structure_cif)

    result, num_chains, num_res = parse_structure(input_structure_cif)

    assert result == structure
    assert num_chains == 2
    assert num_res == 116
