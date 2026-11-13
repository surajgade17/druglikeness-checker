"""PAINS and Brenk structural alert matching."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import FilterCatalog


def find_alerts(molecule: Chem.Mol) -> list[dict[str, str]]:
    """Return matched PAINS and Brenk alert names and available SMARTS patterns."""
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")
    found: list[dict[str, str]] = []
    for catalog_kind, label in ((FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A, "PAINS"),
                                (FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK, "Brenk")):
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(catalog_kind)
        catalog = FilterCatalog.FilterCatalog(params)
        for entry in catalog.GetMatches(molecule):
            try:
                pattern = entry.GetProp("Smarts") if entry.HasProp("Smarts") else "Not provided by RDKit"
            except (AttributeError, RuntimeError):
                pattern = "Not provided by RDKit"
            found.append({"category": label, "alert": entry.GetDescription(), "matched_substructure": pattern})
    return found
