from utils.validation import parse_smiles


def test_valid_smiles_parses() -> None:
    molecule, error = parse_smiles("CCO")
    assert molecule is not None and error is None


def test_invalid_smiles_returns_friendly_error() -> None:
    molecule, error = parse_smiles("not a smiles")
    assert molecule is None and error is not None


def test_empty_smiles_returns_error() -> None:
    molecule, error = parse_smiles("   ")
    assert molecule is None and "non-empty" in error
