"""Calculate common molecular descriptors and optional QED/SA scores."""
from __future__ import annotations

from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, QED


def calculate_properties(molecule: Chem.Mol) -> dict[str, Any]:
    """Return 2D descriptors for a validated RDKit molecule.

    SA score is obtained from RDKit's optional Contrib module when installed. It is
    explicitly marked unavailable otherwise; descriptor calculation still succeeds.
    """
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")
    values: dict[str, Any] = {
        "Molecular weight": round(Descriptors.MolWt(molecule), 2),
        "cLogP": round(Descriptors.MolLogP(molecule), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(molecule), 2),
        "HBD": int(Lipinski.NumHDonors(molecule)),
        "HBA": int(Lipinski.NumHAcceptors(molecule)),
        "Rotatable bonds": int(Lipinski.NumRotatableBonds(molecule)),
        "Rings": int(rdMolDescriptors.Calc NumRings(molecule)),
        "Heavy atoms": int(molecule.GetNumHeavyAtoms()),
        "QED": round(QED.qed(molecule), 3),
    }
    try:
        from rdkit.Contrib.SA_Score import sascorer  # type: ignore[import-not-found]
        values["SA score"] = round(float(sascorer.calculateScore(molecule)), 2)
        values["SA score status"] = "available"
    except (ImportError, AttributeError, OSError, ValueError, RuntimeError):
        values["SA score"] = None
        values["SA score status"] = "unavailable (RDKit Contrib SA_Score not installed)"
    return values
