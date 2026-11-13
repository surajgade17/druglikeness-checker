from __future__ import annotations

from typing import Any

import pandas as pd
from pandas import Series

from core.properties import calculate_properties
from utils.validation import parse_smiles

FREE_TIER_LIMIT = 20


def read_smiles_csv(uploaded_file: Any) -> pd.DataFrame:
    """Read a CSV file that must contain a 'smiles' column."""
    if uploaded_file is None:
        raise ValueError("A CSV file is required.")
    dataframe = pd.read_csv(uploaded_file)
    if "smiles" not in dataframe.columns:
        raise ValueError("Uploaded CSV must contain a 'smiles' column.")
    return dataframe[["smiles"]].copy()


def process_batch_smiles(smiles_series: Series, limit: int = FREE_TIER_LIMIT) -> pd.DataFrame:
    """Process a series of SMILES strings and return structured results.

    The free tier is explicitly limited to 20 molecules; the app should surface a
    friendly message when the user exceeds the limit.
    """
    if smiles_series is None:
        raise ValueError("A non-empty series of SMILES strings is required.")

    if len(smiles_series) > limit:
        raise ValueError(f"Free tier supports up to {limit} molecules. Please reduce the dataset or use Pro.")

    rows: list[dict[str, Any]] = []
    for index, value in smiles_series.items():
        smiles = str(value).strip()
        molecule, error = parse_smiles(smiles)
        if molecule is None or error:
            rows.append(
                {
                    "index": index,
                    "smiles": smiles,
                    "status": "invalid",
                    "error": error or "Invalid SMILES",
                }
            )
            continue

        properties = calculate_properties(molecule)
        rows.append(
            {
                "index": index,
                "smiles": smiles,
                "status": "ok",
                "Molecular weight": properties.get("Molecular weight"),
                "cLogP": properties.get("cLogP"),
                "TPSA": properties.get("TPSA"),
                "HBD": properties.get("HBD"),
                "HBA": properties.get("HBA"),
                "Rotatable bonds": properties.get("Rotatable bonds"),
                "Rings": properties.get("Rings"),
                "Heavy atoms": properties.get("Heavy atoms"),
                "QED": properties.get("QED"),
                "SA score": properties.get("SA score"),
                "error": "",
            }
        )

    return pd.DataFrame(rows)
