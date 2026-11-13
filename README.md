import pytest
from rdkit import Chem

from core.rules import egan, evaluate_rules, ghose, lipinski, muegge, veber


@pytest.mark.parametrize(
    ("mw, logp, hbd, hba, expected"),
    [
        (180.0, 1.5, 1, 4, True),
        (600.0, 6.0, 6, 12, False),
    ],
)
def test_lipinski(mw: float, logp: float, hbd: int, hba: int, expected: bool) -> None:
    assert lipinski(mw, logp, hbd, hba).passed is expected


@pytest.mark.parametrize(
    ("rotatable_bonds, tpsa, expected"),
    [
        (2, 50.0, True),
        (11, 150.0, False),
    ],
)
def test_veber(rotatable_bonds: int, tpsa: float, expected: bool) -> None:
    assert veber(rotatable_bonds, tpsa).passed is expected


@pytest.mark.parametrize(
    ("mw, logp, mr, atoms, expected"),
    [
        (300.0, 2.0, 80.0, 25, True),
        (100.0, 8.0, 20.0, 10, False),
    ],
)
def test_ghose(mw: float, logp: float, mr: float, atoms: int, expected: bool) -> None:
    assert ghose(mw, logp, mr, atoms).passed is expected


@pytest.mark.parametrize(
    ("logp, tpsa, expected"),
    [
        (2.0, 100.0, True),
        (6.0, 140.0, False),
    ],
)
def test_egan(logp: float, tpsa: float, expected: bool) -> None:
    assert egan(logp, tpsa).passed is expected


@pytest.mark.parametrize(
    ("mw, logp, tpsa, hbd, hba, rotors, rings, expected"),
    [
        (350.0, 2.0, 90.0, 2, 5, 6, 3, True),
        (700.0, 7.0, 180.0, 7, 12, 16, 8, False),
    ],
)
def test_muegge(
    mw: float,
    logp: float,
    tpsa: float,
    hbd: int,
    hba: int,
    rotors: int,
    rings: int,
    expected: bool,
) -> None:
    assert muegge(mw, logp, tpsa, hbd, hba, rotors, rings).passed is expected


def test_evaluate_rules_aspirin() -> None:
    molecule = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert molecule is not None
    results = evaluate_rules(molecule)
    assert [result.name for result in results] == ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]
    assert results[0].passed
    assert results[1].passed
    assert results[3].passed


def test_rules_reject_missing_molecule() -> None:
    with pytest.raises(ValueError):
        evaluate_rules(None)
