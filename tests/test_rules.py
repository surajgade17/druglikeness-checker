import pytest
from rdkit import Chem

from core.rules import egan, evaluate_rules, ghose, lipinski, muegge, veber


@pytest.mark.parametrize("mw,logp,hbd,hba,expected", [(180.0, 1.5, 1, 4, True), (600.0, 6.0, 6, 12, False)])
def test_lipinski(mw: float, logp: float, hbd: int, hba: int, expected: bool) -> None:
    assert lipinski(mw, logp, hbd, hba).passed is expected


@pytest.mark.parametrize("rotors,tpsa,expected", [(2, 50.0, True), (11, 150.0, False)])
def test_veber(rotors: int, tpsa: float, expected: bool) -> None:
    assert veber(rotors, tpsa).passed is expected


@pytest.mark.parametrize("mw,logp,mr,atoms,expected", [(300.0, 2.0, 80.0, 25, True), (100.0, 8.0, 20.0, 10, False)])
def test_ghose(mw: float, logp: float, mr: float, atoms: int, expected: bool) -> None:
    assert ghose(mw, logp, mr, atoms).passed is expected


@pytest.mark.parametrize("logp,tpsa,expected", [(2.0, 100.0, True), (6.0, 140.0, False)])
def test_egan(logp: float, tpsa: float, expected: bool) -> None:
    assert egan(logp, tpsa).passed is expected


@pytest.mark.parametrize("mw,logp,tpsa,hbd,hba,rotors,rings,expected", [(350.0, 2.0, 90.0, 2, 5, 6, 3, True), (700.0, 7.0, 180.0, 7, 12, 16, 8, False)])
def test_muegge(mw: float, logp: float, tpsa: float, hbd: int, hba: int, rotors: int, rings: int, expected: bool) -> None:
    assert muegge(mw, logp, tpsa, hbd, hba, rotors, rings).passed is expected


def test_all_rules_run_on_aspirin_gefitinib_caffeine() -> None:
    examples = ["CC(=O)OC1=CC=CC=C1C(=O)O", "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6", "Cn1c(=O)n(C)c(=O)c2ncn(C)c12"]
    for smiles in examples:
        molecule = Chem.MolFromSmiles(smiles)
        assert molecule is not None
        results = evaluate_rules(molecule)
        assert [item.name for item in results] == ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]


def test_rules_reject_missing_molecule() -> None:
    with pytest.raises(ValueError):
        evaluate_rules(None)
