"""Focused tests for each drug-likeness rule."""
import pytest
from rdkit import Chem

from core.rules import egan, evaluate_rules, ghose, lipinski, muegge, veber


@pytest.mark.parametrize("args, expected", [((180, 2, 1, 4), True), ((600, 6, 6, 12), False)])
def test_lipinski(args, expected):
    assert lipinski(*args).passed is expected

@pytest.mark.parametrize("args, expected", [((2, 50), True), ((11, 150), False)])
def test_veber(args, expected):
    assert veber(*args).passed is expected

@pytest.mark.parametrize("args, expected", [((300, 2, 80, 40), True), ((100, 8, 20, 10), False)])
def test_ghose(args, expected):
    assert ghose(*args).passed is expected

@pytest.mark.parametrize("args, expected", [((2, 100), True), ((6, 140), False)])
def test_egan(args, expected):
    assert egan(*args).passed is expected

@pytest.mark.parametrize("args, expected", [((350, 2, 90, 2, 5, 6, 3), True), ((700, 7, 180, 7, 12, 16, 8), False)])
def test_muegge(args, expected):
    assert muegge(*args).passed is expected


def test_evaluate_rules_aspirin():
    molecule = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    results = evaluate_rules(molecule)
    assert [result.name for result in results] == ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]
    assert results[0].passed and results[1].passed and results[3].passed
    assert not results[2].passed and not results[4].passed


def test_rules_reject_missing_molecule():
    with pytest.raises(ValueError):
        evaluate_rules(None)
