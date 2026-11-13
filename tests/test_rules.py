from __future__ import annotations

from typing import Optional

from rdkit import Chem


def parse_smiles(smiles: object) -> tuple[Optional[Chem.Mol], Optional[str]]:
    """Parse a SMILES string and return a friendly validation error when needed."""
    if not isinstance(smiles, str):
        return None, "Enter a valid SMILES string."

    cleaned = smiles.strip()
    if not cleaned:
        return None, "Enter a non-empty SMILES string."

    try:
        molecule = Chem.MolFromSmiles(cleaned)
    except (TypeError, ValueError, RuntimeError) as exc:
        return None, f"This SMILES could not be parsed: {exc}"

    if molecule is None:
        return None, "This does not appear to be a valid SMILES string. Check the notation and try again."

    return molecule, None
