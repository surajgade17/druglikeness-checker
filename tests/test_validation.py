import pandas as pd
import pytest

from core.alerts import find_alerts
from core.batch import process_batch_smiles
from utils.validation import parse_smiles


@pytest.mark.parametrize(
    "smiles",
    [
        "CCO",
        "Cn1c(=O)n(C)c(=O)c2ncn(C)c12",
        "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6",
    ],
)
def test_parse_smiles_accepts_valid_examples(smiles: str) -> None:
    molecule, error = parse_smiles(smiles)
    assert molecule is not None
    assert error is None


def test_parse_smiles_rejects_empty_and_garbage() -> None:
    for value in ("", "   ", "not a smiles", "x" * 100_001):
        molecule, error = parse_smiles(value)
        assert molecule is None
        assert error is not None


def test_alert_detection_for_known_patterns() -> None:
    nitrobenzene, _ = parse_smiles("[N+](=O)([O-])c1ccccc1")
    assert nitrobenzene is not None
    assert any(alert["category"] == "PAINS" for alert in find_alerts(nitrobenzene))

    michael_acceptor, _ = parse_smiles("C=CC(=O)C")
    assert michael_acceptor is not None
    assert any(alert["category"] == "Brenk" for alert in find_alerts(michael_acceptor))


def test_process_batch_smiles_handles_valid_and_invalid_rows() -> None:
    results = process_batch_smiles(pd.Series(["CCO", "invalid" , "Cn1c(=O)n(C)c(=O)c2ncn(C)c12"]))
    assert len(results) == 3
    assert results["status"].tolist() == ["ok", "invalid", "ok"]
    assert "error" in results.columns


def test_process_batch_smiles_rejects_over_limit() -> None:
    with pytest.raises(ValueError):
        process_batch_smiles(pd.Series(["CCO"] * 21), limit=20)
