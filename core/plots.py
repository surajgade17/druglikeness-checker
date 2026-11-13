from __future__ import annotations

from typing import Any

from rdkit import Chem


ALERT_PATTERNS: list[dict[str, str]] = [
    {
        "category": "PAINS",
        "alert": "Nitroaromatic",
        "smarts": "c1ccccc1[N+](=O)[O-]",
    },
    {
        "category": "PAINS",
        "alert": "Azo",
        "smarts": "[#6]-[#7]=[#7]-[#6]",
    },
    {
        "category": "PAINS",
        "alert": "Acyl halide",
        "smarts": "[#6](=[#8])-[#9,#17,#35,#53]",
    },
    {
        "category": "Brenk",
        "alert": "Michael acceptor",
        "smarts": "[#6]=[#6]-[#6](=[#8])-[#6]",
    },
    {
        "category": "Brenk",
        "alert": "Nitro group",
        "smarts": "[$([N+](=O)[O-])]",
    },
    {
        "category": "Brenk",
        "alert": "Aromatic nitroso",
        "smarts": "c[N+](=O)[O-]",
    },
]


def find_alerts(molecule: Chem.Mol) -> list[dict[str, str]]:
    """Return a list of matched structural alerts with the SMARTS pattern used."""
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")

    matches: list[dict[str, str]] = []
    for entry in ALERT_PATTERNS:
        pattern = Chem.MolFromSmarts(entry["smarts"])
        if pattern is None:
            continue
        if molecule.HasSubstructMatch(pattern):
            matches.append(
                {
                    "category": entry["category"],
                    "alert": entry["alert"],
                    "matched_substructure": entry["smarts"],
                }
            )
    return matches
