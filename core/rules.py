from __future__ import annotations

from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors


def calculate_properties(molecule: Chem.Mol) -> dict[str, Any]:
    """Return common 2D molecular descriptors used in the app."""
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")

    values: dict[str, Any] = {
        "Molecular weight": round(float(Descriptors.MolWt(molecule)), 2),
        "cLogP": round(float(Descriptors.MolLogP(molecule)), 2),
        "TPSA": round(float(rdMolDescriptors.CalcTPSA(molecule)), 2),
        "HBD": int(Lipinski.NumHDonors(molecule)),
        "HBA": int(Lipinski.NumHAcceptors(molecule)),
        "Rotatable bonds": int(Lipinski.NumRotatableBonds(molecule)),
        "Rings": int(rdMolDescriptors.CalcNumRings(molecule)),
        "Heavy atoms": int(molecule.GetNumHeavyAtoms()),
        "QED": round(float(QED.qed(molecule)), 3),
    }

    try:
        from rdkit.Contrib.SA_Score import sascorer  # type: ignore[import-not-found]

        values["SA score"] = round(float(sascorer.calculateScore(molecule)), 2)
        values["SA score status"] = "available"
    except (ImportError, AttributeError, OSError, ValueError, RuntimeError):
        values["SA score"] = None
        values["SA score status"] = "unavailable"

    return values
