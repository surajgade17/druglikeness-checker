from __future__ import annotations

from typing import Any

from rdkit import Chem


def parse_smiles(smiles: Any) -> tuple[Chem.Mol | None, str | None]:
    """Safely parse a SMILES string and return a friendly error if invalid."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None, "SMILES cannot be empty."

    try:
        molecule = Chem.MolFromSmiles(smiles.strip())
    except (TypeError, ValueError, RuntimeError):
        molecule = None

    if molecule is None:
        return None, "That SMILES could not be parsed. Check the notation and try again."
    return molecule, None
