"""DrugLikeness Checker: simple Streamlit UI."""
from __future__ import annotations

from io import BytesIO
from typing import Any

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from core.alerts import find_alerts
from core.batch import FREE_TIER_LIMIT, process_batch_smiles, read_smiles_csv
from core.plots import build_bar_chart, build_radar_chart
from core.properties import calculate_properties
from core.rules import evaluate_rules
from utils.validation import parse_smiles

EXAMPLE_SMILES: dict[str, str] = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Gefitinib": "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6",
    "Erlotinib": "C#CC1=CC=CC(=C1)NC1=NC=NC2=C1C=CC(=C2)OCCOC",
    "Caffeine": "Cn1c(=O)n(C)c(=O)c2ncn(C)c12",
}


@st.cache_data(show_spinner=False)
def analyze_smiles(smiles: str) -> dict[str, Any]:
    """Run analysis for a single SMILES input."""
    molecule, error = parse_smiles(smiles)
    if molecule is None or error:
        return {"error": error or "Invalid SMILES."}

    properties = calculate_properties(molecule)
    rules = evaluate_rules(molecule)
    alerts = find_alerts(molecule)
    return {
        "molecule": molecule,
        "properties": properties,
        "rules": [
            {"name": result.name, "passed": result.passed, "details": result.details, "violations": result.violations}
            for result in rules
        ],
        "alerts": alerts,
        "canonical_smiles": Chem.MolToSmiles(molecule),
    }


@st.cache_data(show_spinner=False)
def prepare_molecule_image(smiles: str) -> bytes:
    """Render a PNG image of a molecule."""
    molecule, error = parse_smiles(smiles)
    if molecule is None:
        raise ValueError(error or "Invalid SMILES provided to image generation.")
    png_image = Draw.MolToImage(molecule, size=(500, 350))
    buffer = BytesIO()
    png_image.save(buffer, format="PNG")
    return buffer.getvalue()


def render_single_molecule() -> None:
    """Single-molecule analysis UI."""
    st.subheader("Single molecule")
    example_names = ["Custom input", *list(EXAMPLE_SMILES.keys())]
    selected = st.selectbox("Use an example molecule", example_names, index=0)
    default_smiles = EXAMPLE_SMILES.get(selected, "") if selected != "Custom input" else "CCO"
    smiles_input = st.text_area("SMILES", value=default_smiles, height=120, help="Paste a valid SMILES string.")

    if st.button("Analyze molecule", type="primary") or smiles_input.strip():
        result = analyze_smiles(smiles_input)
        if "error" in result:
            st.error(result["error"])
            return

        image_bytes = prepare_molecule_image(smiles_input)
        st.image(image_bytes, caption="2D structure", use_container_width=True)
        st.caption(f"Canonical SMILES: `{result['canonical_smiles']}`")

        props = result["properties"]
        rule_results = result["rules"]
        alerts = result["alerts"]

        for index in range(0, len(props.items()), 2):
            cols = st.columns(2)
            items = list(props.items())[index:index + 2]
            for col, (name, value) in zip(cols, items):
                display = value if value is not None else "N/A"
                col.metric(name, f"{display}")

        st.markdown("### Drug-likeness rules")
        for rule in rule_results:
            status = "✅ PASS" if rule["passed"] else "❌ FAIL"
            st.markdown(f"{status} **{rule['name']}** — {rule['details']}")

        st.markdown("### Alerts")
        if alerts:
            for entry in alerts:
                st.warning(f"{entry['category']}: {entry['alert']} | substructure: `{entry['matched_substructure']}`")
        else:
            st.success("No PAINS or Brenk alerts were matched.")

        st.markdown("### Comparison chart")
        radar_chart = build_radar_chart(props)
        bar_chart = build_bar_chart(props)
        st.plotly_chart(radar_chart, use_container_width=True)
        st.plotly_chart(bar_chart, use_container_width=True)


def render_batch_mode() -> None:
    """Batch CSV UI."""
    st.subheader("Batch mode")
    st.caption(f"Upload a CSV with a `smiles` column. Free tier limit: {FREE_TIER_LIMIT} molecules.")
    uploaded = st.file_uploader("Choose CSV", type=["csv"])
    if uploaded is None:
        return

    try:
        raw_df = read_smiles_csv(uploaded)
        results = process_batch_smiles(raw_df["smiles"], limit=FREE_TIER_LIMIT)
    except ValueError as exc:
        st.error(str(exc))
        return

    st.dataframe(results, use_container_width=True)
    csv_data = results.to_csv(index=False).encode("utf-8")
    st.download_button("Download results CSV", csv_data, "druglikeness_results.csv", "text/csv", use_container_width=True)


def main() -> None:
    """Main app entry point."""
    st.title("💊 DrugLikeness Checker")
    st.caption("Quick molecular property screening from SMILES")
    st.info("For research and educational use only. Not medical advice.")
    render_single_molecule()
    render_batch_mode()


if __name__ == "__main__":
    main()
