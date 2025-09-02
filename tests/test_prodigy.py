import json
import tarfile
import tempfile
from io import BufferedReader, TextIOWrapper
from os.path import basename, splitext
from pathlib import Path

import pytest
from Bio.PDB.Model import Model
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from prodigy_prot.modules.parsers import validate_structure
from prodigy_prot.modules.prodigy import (
    Prodigy,
    analyse_contacts,
    analyse_nis,
    calculate_ic,
)

from . import TEST_DATA


@pytest.fixture
def input_model():
    input_f = Path(TEST_DATA, "2oob.pdb")
    parser = PDBParser()
    structure = parser.get_structure(input_f.stem, input_f)
    assert isinstance(structure, Structure)
    return structure.child_list[0]


@pytest.fixture
def compressed_dataset_f():
    return Path(TEST_DATA, "dataset.tgz")


@pytest.fixture
def expected_dataset_json():
    return Path(TEST_DATA, "dataset.json")


@pytest.fixture
def prodigy_class(input_model):
    yield Prodigy(input_model)


def test_calculate_ic(input_model):

    result = calculate_ic(model=input_model, d_cutoff=5.5)

    assert len(result) == 78

    first_hit: tuple[Residue, Residue] = result[0]

    assert first_hit[0].get_resname() == "ASN"
    assert first_hit[1].get_resname() == "LYS"


def test_calculate_ic_with_selection(input_model):

    result = calculate_ic(model=input_model, d_cutoff=5.5, selection={"A": 0, "B": 1})

    assert len(result) == 78

    first_hit: tuple[Residue, Residue] = result[0]

    assert first_hit[0].get_resname() == "ASN"
    assert first_hit[1].get_resname() == "LYS"


def test_analyse_contacts(input_model):

    res_a = input_model["A"][(" ", 931, " ")]
    res_b = input_model["B"][(" ", 6, " ")]
    contact = (res_a, res_b)

    test_input = [contact]

    result = analyse_contacts(test_input)

    expected_output = {
        "AA": 0.0,
        "PP": 0.0,
        "CC": 0.0,
        "AP": 0.0,
        "CP": 1.0,
        "AC": 0.0,
    }

    assert result == expected_output


def test_analyse_nis():

    test_input = {("B", "ARG", "72"): 0.9}
    apolar, polar, charged = analyse_nis(test_input)

    assert apolar == 0.0
    assert polar == 100.0
    assert charged == 0.0


def test_prodigy_predict(prodigy_class):

    prodigy_class.predict()

    assert prodigy_class.nis_a == pytest.approx(35.5, abs=1.0)
    assert prodigy_class.nis_c == pytest.approx(38.0, abs=1.0)
    assert prodigy_class.ba_val == pytest.approx(-6.2, abs=1.0)

    # This is the actual prediction
    assert prodigy_class.kd_val == pytest.approx(2.7e-5, abs=1e-6)


def test_prodigy_as_dict(prodigy_class):

    result = prodigy_class.as_dict()

    assert isinstance(result, dict)
    assert len(result) == 14


def test_prodigy_print_prediction(prodigy_class):

    outfile = tempfile.NamedTemporaryFile(delete=False)
    assert Path(outfile.name).stat().st_size == 0

    prodigy_class.print_prediction(outfile.name)
    assert Path(outfile.name).stat().st_size != 0

    Path(outfile.name).unlink()


def test_prodigy_print_prediction_quiet(prodigy_class):

    outfile = tempfile.NamedTemporaryFile(delete=False)
    assert Path(outfile.name).stat().st_size == 0

    prodigy_class.print_prediction(outfile.name, True)
    assert Path(outfile.name).stat().st_size != 0

    Path(outfile.name).unlink()


def test_prodigy_print_contacts(input_model, prodigy_class):

    res_a = input_model["A"][(" ", 931, " ")]
    res_b = input_model["B"][(" ", 6, " ")]
    prodigy_class.ic_network = [(res_a, res_b)]

    outfile = tempfile.NamedTemporaryFile(delete=False)
    assert Path(outfile.name).stat().st_size == 0

    prodigy_class.print_contacts(outfile.name)
    assert Path(outfile.name).stat().st_size != 0

    Path(outfile.name).unlink()


def test_print_pymol_script(input_model, prodigy_class):
    res_a = input_model["A"][(" ", 931, " ")]
    res_b = input_model["B"][(" ", 6, " ")]
    prodigy_class.ic_network = [(res_a, res_b)]

    outfile = tempfile.NamedTemporaryFile(delete=False)
    assert Path(outfile.name).stat().st_size == 0

    prodigy_class.print_pymol_script(outfile.name)
    assert Path(outfile.name).stat().st_size != 0

    Path(outfile.name).unlink()


@pytest.mark.integration
def test_dataset_prediction(compressed_dataset_f, expected_dataset_json):
    """
    Test method to compare prediction for 80 dataset cases with
        expected values.
    """
    # load expected data from json
    with open(expected_dataset_json) as fh:
        expected_data = json.load(fh)

    # load dataset PDBs
    dataset = tarfile.open(compressed_dataset_f)
    parser = PDBParser(QUIET=True)

    keys_equal = ["AA", "PP", "CC", "AP", "CP", "AC"]
    diffs = {"ba_val": [], "nis_a": [], "nis_c": []}

    # run prodigy for each dataset in the PDB
    for entry in dataset:
        s_name, s_ext = splitext(basename(entry.name))

        # skip system files in archive
        if not s_name.isalnum() or s_ext != ".pdb":
            continue

        handle = dataset.extractfile(entry)

        # Wrap filehandle to ensure string file handle in Python 3
        handle = TextIOWrapper(BufferedReader(handle))  # type: ignore

        parsed_structure = parser.get_structure(s_name, handle)
        assert isinstance(parsed_structure, Structure)

        models = validate_structure(parsed_structure, selection=["A", "B"])

        # Test for structure object
        # Check if it's a list and all elements are Model objects
        assert isinstance(models, list) and all(
            isinstance(item, Model) for item in models
        )
        # assert isinstance(s, list[Model])

        #  run prediction and retrieve result dict
        for m in models:
            prod = Prodigy(m, selection=["A", "B"])
            prod.predict()
            results = prod.as_dict()

            # check for equality of prdicted interface residues
            for k in keys_equal:
                observed_value = results[k]
                expected_value = expected_data[s_name][k]
                assert observed_value == pytest.approx(expected_value)

            # check that NIS and binding afinity values are within 2% of
            #  expected values and add diffs for summary
            for k in diffs.keys():
                delta = abs(results[k] / expected_data[s_name][k] - 1)
                # assume a difference of less then 2%
                assert delta == pytest.approx(0, abs=0.02)
                diffs[k].append(delta)
