import math
import tempfile
from pathlib import Path

from prodigy_prot.modules.utils import check_path, dg_to_kd


def test_check_path():

    temp_f = tempfile.NamedTemporaryFile(delete=False)

    result = check_path(temp_f.name)

    assert result == temp_f.name

    Path(temp_f.name).unlink()


def test_dg_to_kd():

    assert math.isclose(dg_to_kd(0.0), 1.0, rel_tol=1e-9)
