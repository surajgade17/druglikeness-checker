from __future__ import annotations

from typing import Any

import pandas as pd
from pandas import Series

from core.alerts import find_alerts
from core.properties import calculate_properties
from core.rules import evaluate_rules
from utils.validation import parse_smiles

FREE_TIER_LIMIT = 20


def read_smiles_csv(uploaded_file: Any) -> pd.DataFrame:
    """Read a CSV file and validate that it contains a `smiles` column."""
    if uploaded_file is None:
        raise ValueError("A CSV file is required.")
    try:
        dataframe = pd.read_csv(uploaded_file)
    except (pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise ValueError(f"Could not read CSV: {exc}") from exc

    if "smiles" not in dataframe.columns:
        raise ValueError("Uploaded CSV must contain a 'smiles' column.")
    return dataframe[["smiles"]].copy()


def process_batch_smiles(smiles_series: Series, limit: int = FREE_TIER_LIMIT) -> pd.DataFrame:
    """Process a series of SMILES strings and return results."""
    if smiles_series is None:
        raise ValueError("A non-empty series of SMILES strings is required.")
    if len(smiles_series) > limit:
        raise ValueError(f"Free tier supports up to {limit} molecules. Reduce the dataset and retry.")

    rows: list[dict[str, Any]] = []
    for index, value in smiles_series.items():
        smiles = "" if pd.isna(value) else str(value).strip()
        molecule, error = parse_smiles(smiles)
        if molecule is None:
            rows.append({"index": index, "smiles": smiles, "status": "invalid", "error": error or "Invalid SMILES"})
            continue

        properties = calculate_properties(molecule)
        rules = evaluate_rules(molecule)
        alerts = find_alerts(molecule)
        rows.append(
            {
                "index": index,
                "smiles": smiles,
                "status": "ok",
                **{k: v for k, v in properties.items() if not k.endswith("status")},
                "rules_passed": sum(rule.passed for rule in rules),
                "alert_count": len(alerts),
                "error": "",
            }
        )

    return pd.DataFrame(rows)
