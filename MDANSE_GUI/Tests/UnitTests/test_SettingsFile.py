import pytest

from MDANSE_GUI.Session.StructuredSession import SettingsFile, SettingsGroup


@pytest.fixture(scope="module")
def settings_file(tmp_path):
    temp_name = tmp_path / "settings"
    sf = SettingsFile(temp_name.stem, tmp_path)
    yield sf

def test_loading(settings_file: "SettingsFile"):
    new_group = settings_file.group("blam")
    new_group.add("blob", "000", "useless test variable")
    settings_file.save_values()
    settings_file.load_from_file()
    loaded_group = settings_file.group("blam")
    loaded_value = loaded_group.get("blob")
    assert loaded_value == "000"


def test_reloading(tmp_path):
    temp_name = tmp_path / "output"

    settings_file = SettingsFile(temp_name.stem, tmp_path)
    new_group = settings_file.group("blam")
    new_group.add("blob", "000", "useless test variable")
    settings_file.save_values()

    settings_file2 = SettingsFile(temp_name.stem, tmp_path)
    settings_file2.load_from_file()
    loaded_group = settings_file2.group("blam")
    loaded_value = loaded_group.get("blob")

    assert loaded_value == "000"
