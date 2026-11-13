from __future__ import annotations

import os
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

st.set_page_config(page_title="DrugLikeness Checker", page_icon="💊", layout="wide")

EXAMPLE_SMILES: dict[str, str] = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Gefitinib": "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6",
    "Erlotinib": "COCCOC1=CC=C(C=C1)C2=NC(=NC=C2)NC3=CC=C(C=C3)OCCOC",
}


@st.cache_data
def analyze_smiles(smiles: str) -> dict[str, Any]:
    """Run the chemistry analysis pipeline for a single SMILES string."""
    molecule, error = parse_smiles(smiles)
    if error or molecule is None:
        return {"error": error or "Invalid SMILES."}

    properties = calculate_properties(molecule)
    rules = evaluate_rules(molecule)
    alerts = find_alerts(molecule)

    return {
        "molecule": molecule,
        "properties": properties,
        "rules": [
            {
                "name": result.name,
                "passed": result.passed,
                "details": result.details,
                "violations": result.violations,
            }
            for result in rules
        ],
        "alerts": alerts,
        "formatted_smiles": smiles.strip(),
    }


@st.cache_data
def prepare_molecule_image(smiles: str) -> bytes:
    """Render a 2D structure image to PNG bytes for the given SMILES."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES provided to image generation.")
    png_bytes = Draw.MolToImage(molecule, size=(500, 350))
    buffer = BytesIO()
    png_bytes.save(buffer, format="PNG")
    return buffer.getvalue()


def render_badge(passed: bool) -> str:
    """Return a small HTML badge string for pass/fail display."""
    color = "#22c55e" if passed else "#ef4444"
    label = "PASS" if passed else "FAIL"
    return f'<span style="display:inline-block;padding:0.25rem 0.6rem;border-radius:9999px;background:{color};color:white;font-weight:700;font-size:0.75rem;">{label}</span>'


def render_single_molecule() -> None:
    """Render the single-molecule analysis UI."""
    st.subheader("Single molecule")

    example_key = st.selectbox("Use an example molecule", ["Custom input", *list(EXAMPLE_SMILES.keys())], index=0)
    if example_key != "Custom input":
        default_smiles = EXAMPLE_SMILES[example_key]
    else:
        default_smiles = "CCO"

    smiles_input = st.text_area("SMILES", value=default_smiles, height=120, help="Paste a valid SMILES string.")
    if st.button("Analyze molecule", type="primary") or smiles_input.strip():
        try:
            result = analyze_smiles(smiles_input)
        except ValueError as exc:
            st.error(f"Unexpected error: {exc}")
            return

        if "error" in result:
            st.error(result["error"])
            return

        st.markdown("### Molecule overview")
        molecule = Chem.MolFromSmiles(smiles_input)
        if molecule is not None:
            image_bytes = prepare_molecule_image(smiles_input)
            st.image(image_bytes, caption="2D structure", width=500)

        properties = result["properties"]
        rule_results = result["rules"]
        alerts = result["alerts"]

        property_rows = [
            ("Molecular weight", properties.get("Molecular weight"), "Da"),
            ("cLogP", properties.get("cLogP"), ""),
            ("TPSA", properties.get("TPSA"), "Å²"),
            ("HBD", properties.get("HBD"), ""),
            ("HBA", properties.get("HBA"), ""),
            ("Rotatable bonds", properties.get("Rotatable bonds"), ""),
            ("Rings", properties.get("Rings"), ""),
            ("Heavy atoms", properties.get("Heavy atoms"), ""),
            ("QED", properties.get("QED"), ""),
            ("SA score", properties.get("SA score"), ""),
        ]

        cols = st.columns(2)
        for idx, (label, value, unit) in enumerate(property_rows):
            col = cols[idx % 2]
            if value is None:
                display = "N/A"
            elif isinstance(value, (int, float)) and label != "SA score":
                display = f"{value:g}{unit}"
            else:
                display = f"{value}{unit}" if unit else str(value)
            col.metric(label, display)

        st.markdown("### Drug-likeness rules")
        for rule in rule_results:
            st.markdown(
                f"<div style='margin-bottom:0.75rem;'><strong>{rule['name']}</strong> {render_badge(rule['passed'])} - {rule['details']}</div>",
                unsafe_allow_html=True,
            )

        st.markdown("### Alerts")
        if alerts:
            for alert in alerts:
                st.markdown(
                    f"- <strong>{alert['category']}</strong>: {alert['alert']}  |  substructure: `{alert['matched_substructure']}`",
                    unsafe_allow_html=True,
                )
        else:
            st.success("No PAINS or Brenk alerts were matched.")

        radar_chart = build_radar_chart(properties)
        bar_chart = build_bar_chart(properties)
        st.markdown("### Comparison chart")
        st.plotly_chart(radar_chart, use_container_width=True)
        st.plotly_chart(bar_chart, use_container_width=True)

        st.caption("For research and educational use only. Not medical advice.")


def render_batch_mode() -> None:
    """Render the CSV upload and batch analysis workflow."""
    st.subheader("Batch mode")
    uploaded = st.file_uploader("Upload a CSV with a `smiles` column", type=["csv"])

    if uploaded is None:
        st.info(f"Batch mode is limited to {FREE_TIER_LIMIT} molecules in the free tier.")
        return

    try:
        raw_df = read_smiles_csv(uploaded)
    except ValueError as exc:
        st.error(str(exc))
        return

    try:
        results = process_batch_smiles(raw_df["smiles"], limit=FREE_TIER_LIMIT)
    except ValueError as exc:
        st.error(str(exc))
        return

    st.dataframe(results, use_container_width=True)

    csv_data = results.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download CSV results",
        data=csv_data,
        file_name="druglikeness_batch_results.csv",
        mime="text/csv",
    )


def pro_status() -> bool:
    """Check whether the user has provided a valid Pro key for premium features."""
    configured_key = os.getenv("DRUGLIKENESS_PRO_KEY", "")
    if not configured_key:
        try:
            configured_key = st.secrets["DRUGLIKENESS_PRO_KEY"]  # type: ignore[index]
        except Exception:
            configured_key = ""
    user_key = st.session_state.get("pro_key", "")
    return bool(configured_key and configured_key == user_key)


def render_pro_stub() -> None:
    """Render a gated Pro features panel with a clearly marked TODO for billing."""
    st.subheader("Pro features")
    st.info("TODO: integrate payment flow and licensing backend. This is a stub implementation for gating premium features.")

    entered_key = st.text_input("Pro license key", type="password", key="pro_key")

    if not entered_key:
        st.caption("Use a valid key to enable premium features. Example key format: PRO-XXXX")
        return

    enabled = pro_status()
    if enabled:
        st.success("Pro mode enabled")
        st.markdown(
            "- Up to 5000 molecules per batch\n"
            "- PDF report export (planned)\n"
            "- Tanimoto similarity search against a user-supplied reference set (planned)"
        )
    else:
        st.warning("Invalid or missing license key. Premium features remain locked.")


def main() -> None:
    """Main Streamlit entry point for the app."""
    st.title("DrugLikeness Checker")
    st.markdown("A quick screen for SMILES-based drug-likeness and basic ADMET-style checks.")

    render_single_molecule()
    render_batch_mode()
    render_pro_stub()

    st.caption("For research and educational use only. Not medical advice.")


if __name__ == "__main__":
    main()
