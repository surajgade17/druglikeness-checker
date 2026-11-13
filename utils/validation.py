"""SMILES parsing and input validation utilities."""
from __future__ import annotations

from typing import Optional

from rdkit import Chem


def parse_smiles(smiles: object) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Parse a SMILES string, returning a molecule or a friendly validation error."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None, "Enter a non-empty SMILES string."
    candidate = smiles.strip()
    try:
        molecule = Chem.MolFromSmiles(candidate)
    except (ValueError, RuntimeError) as exc:
        return None, f"This SMILES could not be parsed: {exc}"
    if molecule is None:
        return None, "This does not appear to be a valid SMILES string. Check the notation and try again."
    return molecule, None
