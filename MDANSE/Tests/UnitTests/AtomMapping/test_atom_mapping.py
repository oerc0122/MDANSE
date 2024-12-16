import pytest

from MDANSE.Framework.AtomMapping.atom_mapping import guess_element
from MDANSE.Framework.AtomMapping.atom_mapping import fill_remaining_labels
from MDANSE.Framework.AtomMapping.atom_mapping import get_element_from_mapping
from MDANSE.Framework.AtomMapping.atom_mapping import check_mapping_valid
from MDANSE.Framework.AtomMapping.atom_mapping import AtomLabel


def test_guess_element_hydrogen_to_hydrogen_0():
    h = guess_element("1H")
    assert h == "H"


def test_guess_element_hydrogen_to_hydrogen_1():
    h = guess_element("h")
    assert h == "H"


def test_guess_element_carbon_double_bond_symbol_to_carbon():
    c = guess_element("C=")
    assert c == "C"


def test_guess_element_gold_to_gold():
    au = guess_element("Au")
    assert au == "Au"


def test_guess_element_raise_error_when_unable_to_guess_match():
    with pytest.raises(AttributeError):
        guess_element("aaa")


def test_guess_element_carbon_symbol_to_carbon_12():
    c = guess_element("C", mass=12)
    assert c == "C12"


def test_guess_element_carbon_symbol_to_carbon_13():
    c = guess_element("C", mass=13)
    assert c == "C13"


def test_guess_element_carbon_symbol_to_carbon():
    c = guess_element("C", mass=12.01)
    assert c == "C"


def test_guess_element_carbon_symbol_to_deuterium():
    h2 = guess_element("H", mass=2)
    assert h2 == "H2"


def test_guess_element_sodium_to_sodium():
    na = guess_element("NA")
    assert na == "Na"


def test_guess_element_sodium_with_mass_to_sodium():
    na = guess_element("NA", mass=22.98977)
    assert na == "Na"


def test_guess_element_sod_with_mass_to_sodium():
    na = guess_element("SOD", mass=22.98977)
    assert na == "Na"


def test_guess_element_chlorine_to_chlorine():
    na = guess_element("CL")
    assert na == "Cl"


def test_get_element_from_mapping_with_entry_but_no_kwargs():
    mapping = {"": {"label1": "C"}}
    element = get_element_from_mapping(mapping, "label1")
    assert element == "C"


def test_get_element_from_mapping_with_entry():
    mapping = {"molecule=mol1": {"label1": "C"}}
    element = get_element_from_mapping(mapping, "label1", molecule="mol1")
    assert element == "C"


def test_get_element_from_mapping_with_entry_two_group_labels():
    mapping = {"grp1=1;grp2=1": {"label1": "C"}}
    element = get_element_from_mapping(mapping, "label1", grp1="1", grp2="1")
    assert element == "C"


def test_get_element_from_mapping_with_no_entry():
    element = get_element_from_mapping({}, "C", molecule="mol1")
    assert element == "C"


def test_get_element_from_mapping_with_no_entry_raises_error_with_bad_label():
    with pytest.raises(AttributeError):
        get_element_from_mapping({}, "aaa", molecule="mol1")


def test_fill_remaining_labels_fill_mapping_correctly_already_filled_1():
    labels = [
        AtomLabel("label1", molecule="mol1"),
        AtomLabel("label2", molecule="mol1"),
    ]
    mapping = {"molecule=mol1": {"label1": "C", "label2": "C"}}
    fill_remaining_labels(mapping, labels)
    assert mapping == {"molecule=mol1": {"label1": "C", "label2": "C"}}


def test_fill_remaining_labels_fill_mapping_correctly_already_filled_2():
    labels = [
        AtomLabel("label1", molecule="mol1"),
        AtomLabel("label1", molecule="mol1"),
    ]
    mapping = {"molecule=mol1": {"label1": "C"}}
    fill_remaining_labels(mapping, labels)
    assert mapping == {"molecule=mol1": {"label1": "C"}}


def test_fill_remaining_labels_fill_mapping_correctly_missing_mapping():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1": {"label1": "C"}}
    fill_remaining_labels(mapping, labels)
    assert mapping == {"molecule=mol1": {"label1": "C", "C1": "C"}}


def test_check_mapping_valid_as_elements_are_correct_1():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1": {"label1": "C", "C1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_valid_as_elements_are_correct_2():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1": {"C1": "C", "label1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_valid_as_elements_are_correct_3():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol2")]
    mapping = {"molecule=mol1": {"label1": "C"}, "molecule=mol2": {"C1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_valid_as_elements_are_correct_4():
    labels = [AtomLabel("label1", type="1*"), AtomLabel("C1", type="2*")]
    mapping = {"type=1*": {"label1": "C"}, "type=2*": {"C1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_valid_as_elements_are_correct_5():
    labels = [AtomLabel("label1", mass="12"), AtomLabel("C1", mass="13")]
    mapping = {"mass=12": {"label1": "C"}, "mass=13": {"C1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_valid_as_elements_are_correct_6():
    labels = [AtomLabel("label=1", molecule="mol=1"), AtomLabel("C=1", molecule="mol=2")]
    mapping = {"molecule=mol1": {"label1": "C"}, "molecule=mol2": {"C1": "C"}}
    assert check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_as_elements_not_correct():
    labels = [
        AtomLabel("label1", molecule="mol1"),
        AtomLabel("label2", molecule="mol1"),
    ]
    mapping = {"molecule=mol1": {"label1": "C", "label2": "aaa"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_missing_mappings_1():
    labels = [
        AtomLabel("label1", molecule="mol1"),
        AtomLabel("label2", molecule="mol1"),
    ]
    mapping = {"molecule=mol1": {"label1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_missing_mappings_2():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol2")]
    mapping = {"molecule=mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_incorrect_keys_1():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1": {"label1": "C", "label2": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_incorrect_keys_2():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1": {"label1": "C", "C1": "C", "label2": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_1():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule==mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_2():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule===mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_3():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mol1;": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_4():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {";molecule=mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_5():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule;=mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_6():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"mole;cule=mol1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_7():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol1")]
    mapping = {"molecule=mo;l1": {"label1": "C", "C1": "C"}}
    assert not check_mapping_valid(mapping, labels)


def test_check_mapping_not_valid_due_to_bad_keys_8():
    labels = [AtomLabel("label1", molecule="mol1"), AtomLabel("C1", molecule="mol2")]
    mapping = {"molecule=mol1": {"label1": "C"}, "molec;ule=mol2": {"C1": "C"}}
    assert not check_mapping_valid(mapping, labels)
