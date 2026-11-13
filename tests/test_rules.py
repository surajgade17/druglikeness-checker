from __future__ import annotations

from typing import Any

from rdkit import Chem

MAX_SMILES_LENGTH = 100_000


def parse_smiles(smiles: Any) -> tuple[Chem.Mol | None, str | None]:
    """Safely parse a SMILES string and return a helpful validation error."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None, "SMILES cannot be empty."

    text = smiles.strip()
    if len(text) > MAX_SMILES_LENGTH:
        return None, f"SMILES is too long (maximum {MAX_SMILES_LENGTH} characters)."

    try:
        molecule = Chem.MolFromSmiles(text)
    except (TypeError, ValueError, RuntimeError):
        molecule = None

    if molecule is None:
        return None, "That SMILES could not be parsed. Check the notation and try again."
    return molecule, None
