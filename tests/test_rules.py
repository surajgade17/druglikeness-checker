import pandas as pd

from core.batch import process_batch_smiles
from utils.validation import parse_smiles


def test_parse_smiles_handles_long_garbage() -> None:
    molecule, error = parse_smiles("x" * 4096)
    assert molecule is None and error is not None


def test_process_batch_smiles_validates_expected_input() -> None:
    series = pd.Series(["CCO", "CCN", "not a smiles"])
    result = process_batch_smiles(series, limit=10)
    assert len(result) == 3
    assert result["status"].tolist() == ["ok", "ok", "invalid"]


def test_process_batch_smiles_rejects_items_over_limit() -> None:
    series = pd.Series(["CCO"] * 21)
    try:
        process_batch_smiles(series, limit=20)
        raise AssertionError("Expected ValueError for oversized batch")
    except ValueError:
        pass


def test_parse_smiles_handles_caffeine() -> None:
    molecule, error = parse_smiles("Cn1c(=O)n(C)c(=O)c2ncn(C)c12")
    assert molecule is not None and error is None


def test_parse_smiles_handles_gefitinib() -> None:
    molecule, error = parse_smiles("COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6")
    assert molecule is not None and error is None


def test_parse_smiles_handles_garbage() -> None:
    molecule, error = parse_smiles("smiles,CCO")
    assert molecule is None and error is not None
