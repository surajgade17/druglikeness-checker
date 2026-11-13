import pandas as pd

from core.alerts import find_alerts
from core.batch import process_batch_smiles
from utils.validation import parse_smiles


def test_smiles_valid_invalid_empty_and_long() -> None:
    for value in ("CCO", "Cn1c(=O)n(C)c(=O)c2ncn(C)c12"):
        molecule, error = parse_smiles(value)
        assert molecule is not None and error is None
    for value in ("not a smiles", "", "  ", "x" * 100_001):
        molecule, error = parse_smiles(value)
        assert molecule is None and error is not None


def test_pa_ins_and_brenk_alert_patterns() -> None:
    nitrobenzene, _ = parse_smiles("[N+](=O)([O-])c1ccccc1")
    michael_acceptor, _ = parse_smiles("C=CC(=O)C")
    assert nitrobenzene is not None and any(x["category"] == "PAINS" for x in find_alerts(nitrobenzene))
    assert michael_acceptor is not None and any(x["category"] == "Brenk" for x in find_alerts(michael_acceptor))


def test_batch_results_include_invalid_reason() -> None:
    results = process_batch_smiles(pd.Series(["CCO", "garbage"]))
    assert results["status"].tolist() == ["ok", "invalid"]
    assert results.iloc[1]["error"]
