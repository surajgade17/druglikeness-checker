"""DrugLikeness Checker Streamlit interface."""
from __future__ import annotations

from io import BytesIO
from typing import Any

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
}


@st.cache_data(show_spinner=False)
def analyze_smiles(smiles: str) -> dict[str, Any]:
    """Calculate descriptors, rule outcomes, and structural alerts."""
    molecule, error = parse_smiles(smiles)
    if molecule is None:
        return {"error": error or "Invalid SMILES."}
    rules = evaluate_rules(molecule)
    return {
        "molecule": molecule,
        "properties": calculate_properties(molecule),
        "rules": [
            {"name": item.name, "passed": item.passed, "details": item.details, "violations": item.violations}
            for item in rules
        ],
        "alerts": find_alerts(molecule),
        "canonical_smiles": Chem.MolToSmiles(molecule),
    }


@st.cache_data(show_spinner=False)
def prepare_molecule_image(smiles: str) -> bytes:
    """Return a PNG drawing for valid SMILES."""
    molecule, error = parse_smiles(smiles)
    if molecule is None:
        raise ValueError(error or "Invalid SMILES.")
    image = Draw.MolToImage(molecule, size=(500, 350))
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    return buffer.getvalue()


def render_single_molecule() -> None:
    """Render single-molecule input and results."""
    st.subheader("Single molecule")
    names = ["Custom input", *EXAMPLE_SMILES]
    selected = st.selectbox("Use an example molecule", names)
    default = EXAMPLE_SMILES.get(selected, "") if selected != "Custom input" else "CCO"
    smiles = st.text_area("SMILES", value=default, height=100, help="Enter a valid SMILES string.")
    if not st.button("Analyze molecule", type="primary") and not smiles.strip():
        return

    result = analyze_smiles(smiles)
    if "error" in result:
        st.error(result["error"])
        return

    st.markdown("### Molecule overview")
    st.image(prepare_molecule_image(smiles), caption="2D molecular structure", use_container_width=True)
    st.caption(f"Canonical SMILES: `{result['canonical_smiles']}`")
    properties = result["properties"]
    entries = [(key, value) for key, value in properties.items() if not key.endswith("status")]
    for start in range(0, len(entries), 2):
        cols = st.columns(2)
        for col, (label, value) in zip(cols, entries[start:start + 2]):
            col.metric(label, "N/A" if value is None else value)

    st.markdown("### Drug-likeness rules")
    for rule in result["rules"]:
        st.markdown(f"{'✅ PASS' if rule['passed'] else '❌ FAIL'} **{rule['name']}** — {rule['details']}")

    st.markdown("### Structural alerts")
    if result["alerts"]:
        for alert in result["alerts"]:
            st.warning(f"{alert['category']}: {alert['alert']}")
    else:
        st.success("No PAINS or Brenk alerts were matched.")

    st.markdown("### Comparison charts")
    st.plotly_chart(build_radar_chart(properties), use_container_width=True)
    st.plotly_chart(build_bar_chart(properties), use_container_width=True)


def render_batch_mode() -> None:
    """Render CSV batch analysis and download."""
    st.subheader("Batch mode")
    st.caption(f"Upload a CSV with a `smiles` column. Limit: {FREE_TIER_LIMIT} molecules.")
    uploaded = st.file_uploader("Upload CSV", type=["csv"])
    if uploaded is None:
        return
    try:
        frame = read_smiles_csv(uploaded)
        results = process_batch_smiles(frame["smiles"], limit=FREE_TIER_LIMIT)
    except (ValueError, UnicodeDecodeError) as exc:
        st.error(str(exc))
        return
    st.dataframe(results, use_container_width=True, hide_index=True)
    st.download_button(
        "Download CSV results", results.to_csv(index=False).encode("utf-8"),
        "druglikeness_results.csv", "text/csv",
    )


def main() -> None:
    """Render the application."""
    st.title("💊 DrugLikeness Checker")
    st.caption("A quick screen for SMILES-based molecular properties and drug-likeness filters.")
    st.info("For research and educational use only. Not medical advice.")
    render_single_molecule()
    st.divider()
    render_batch_mode()
    st.caption("Drug-likeness filters are heuristics, not efficacy, toxicity, or safety predictions.")


if __name__ == "__main__":
    main()
