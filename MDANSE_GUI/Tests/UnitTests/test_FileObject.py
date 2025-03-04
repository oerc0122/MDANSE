import pytest

from MDANSE_GUI.DataViewModel.TrajectoryHolder import FileObject

REFERENCE_BYTES = b"TeStCaSeFoRtHeCaChEfuNcTiOn"
REFERENCE_HASH = "7105425fa73f3e6a31b72ae2ee36235bcf1f4883b16e674b80a04d6587d995ec"


@pytest.fixture(scope="module")
def temporary_fileobject(tmp_path):
    fname = tmp_path / "fileobj"
    with open(fname, "wb") as target:
        target.write(REFERENCE_BYTES)
    fob = FileObject()
    fob.setFilename(fname)
    return fob


def test_hash_value(temporary_fileobject):
    file_hash = temporary_fileobject.hash
    assert file_hash == REFERENCE_HASH
