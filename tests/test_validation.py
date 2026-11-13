"""Tests for safe SMILES parsing."""
from utils.validation import parse_smiles


def test_valid_smiles_parses():
    molecule, error = parse_smiles("CCO")
    assert molecule is not None and error is None


def test_invalid_smiles_returns_friendly_error():
    molecule, error = parse_smiles("not a smiles")
    assert molecule is None and error


def test_empty_smiles_returns_error():
    molecule, error = parse_smiles("  ")
    assert molecule is None and "non-empty" in error
