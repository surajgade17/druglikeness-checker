"""DrugLikeness Checker: lightweight Streamlit molecular screening demo."""
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

st.set_page_config(page_title="DrugLikeness Checker", page_icon="💊", layout="centered")

EXAMPLE_SMILES: dict[str, str] = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Gefitinib": "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6",
    "Erlotinib": "C#CC1=CC=CC(=C1)NC1=NC=NC2=C1C=CC(=C2)OCCOC",
    "Caffeine": "Cn1c(=O)n(C)c(=O)c2ncn(C)c12",
}


@st.cache_data(show_spinner=False)
def analyze_smiles(smiles: str) -> dict[str, Any]:
    """Calculate descriptors, rules, and structural alerts for one molecule."""
    molecule, error = parse_smiles(smiles)
    if molecule is None:
        return {"error": error or "Invalid SMILES."}

    return {
        "molecule": molecule,
        "properties": calculate_properties(molecule),
        "rules": [
            {"name": item.name, "passed": item.passed, "details": item.details, "violations": item.violations}
            for item in evaluate_rules(molecule)
        ],
        "alerts": find_alerts(molecule),
        "canonical_smiles": Chem.MolToSmiles(molecule),
    }


@st.cache_data(show_spinner=False)
def prepare_molecule_image(smiles: str) -> bytes:
    """Render a PNG structure image."""
    molecule, error = parse_smiles(smiles)
    if molecule is None:
        raise ValueError(error or "Invalid SMILES.")
    image = Draw.MolToImage(molecule, size=(500, 350))
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    return buffer.getvalue()


def render_single_molecule() -> None:
    """Single-molecule analysis UI."""
    st.subheader("Single molecule")
    example_names = ["Custom input", *EXAMPLE_SMILES.keys()]
    selected = st.selectbox("Try an example molecule", example_names)
    default_smiles = EXAMPLE_SMILES.get(selected, "") if selected != "Custom input" else "CCO"
    smiles = st.text_area("SMILES", value=default_smiles, height=110)

    if not st.button("Analyze molecule", type="primary") and not smiles.strip():
        return

    result = analyze_smiles(smiles)
    if "error" in result:
        st.error(result["error"])
        return

    st.markdown("### Structure")
    st.image(prepare_molecule_image(smiles), caption="2D structure", use_container_width=True)
    st.caption(f"Canonical SMILES: `{result['canonical_smiles']}`")

    properties = result["properties"]
    property_items = list(properties.items())
    for start in range(0, len(property_items), 2):
        cols = st.columns(2)
        for col, (label, value) in zip(cols, property_items[start:start + 2]):
            col.metric(label, "N/A" if value is None else value)

    st.markdown("### Drug-likeness rules")
    for rule in result["rules"]:
        status = "✅ PASS" if rule["passed"] else "❌ FAIL"
        st.markdown(f"{status} **{rule['name']}** — {rule['details']}")

    st.markdown("### Alerts")
    if result["alerts"]:
        for alert in result["alerts"]:
            st.warning(f"{alert['category']}: {alert['alert']}")
    else:
        st.success("No PAINS or Brenk alerts were matched.")

    st.markdown("### Comparison charts")
    st.plotly_chart(build_radar_chart(properties), use_container_width=True)
    st.plotly_chart(build_bar_chart(properties), use_container_width=True)


def render_batch_mode() -> None:
    """Batch CSV input and results UI."""
    st.subheader("Batch CSV")
    st.caption(f"Upload a CSV with a `smiles` column. Free tier limit: {FREE_TIER_LIMIT} molecules.")
    uploaded = st.file_uploader("Choose CSV", type=["csv"])
    if uploaded is None:
        return

    try:
        frame = read_smiles_csv(uploaded)
        results = process_batch_smiles(frame["smiles"], limit=FREE_TIER_LIMIT)
    except (ValueError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        st.error(str(exc))
        return

    st.dataframe(results, use_container_width=True, hide_index=True)
    st.download_button(
        "Download results CSV",
        results.to_csv(index=False).encode("utf-8"),
        "druglikeness_results.csv",
        "text/csv",
        use_container_width=True,
    )


def main() -> None:
    """App entry point."""
    st.title("💊 DrugLikeness Checker")
    st.caption("Quick molecular property screening from SMILES.")
    st.info("For research and educational use only. Not medical advice.")

    render_single_molecule()
    st.divider()
    render_batch_mode()


if __name__ == "__main__":
    main()
