"""DrugLikeness Checker: lightweight Streamlit molecular screening demo."""
from __future__ import annotations

from io import BytesIO
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Lipinski, QED, rdMolDescriptors, FilterCatalog

st.set_page_config(page_title="DrugLikeness Checker", page_icon="💊", layout="centered")
st.title("💊 DrugLikeness Checker")
st.caption("Quick molecular property screening from SMILES")
st.info("For research and educational use only. Not medical advice.")

EXAMPLES = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Gefitinib": "COC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCOCC3)C4=CC5=C(C=C4)N=CC(=C5)C6=CC=CC=C6",
    "Erlotinib": "C#CC1=CC=CC(=C1)NC1=NC=NC2=C1C=CC(=C2)OCCOC",
}

@st.cache_data(show_spinner=False)
def analyze(smiles: str) -> dict[str, Any]:
    """Calculate descriptors and filter results for a SMILES string."""
    if not smiles or not smiles.strip():
        return {"error": "Enter a SMILES string first."}
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
    except (ValueError, RuntimeError, TypeError):
        mol = None
    if mol is None:
        return {"error": "That SMILES could not be parsed. Check the notation and try again."}

    mw = float(Descriptors.MolWt(mol))
    logp = float(Descriptors.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd, hba = int(Lipinski.NumHDonors(mol)), int(Lipinski.NumHAcceptors(mol))
    rotors = int(Lipinski.NumRotatableBonds(mol))
    rings, heavy = int(rdMolDescriptors.CalcNumRings(mol)), int(mol.GetNumHeavyAtoms())
    props = {"Molecular weight (Da)": round(mw, 2), "cLogP": round(logp, 2), "TPSA (Å²)": round(tpsa, 2),
             "HBD": hbd, "HBA": hba, "Rotatable bonds": rotors, "Rings": rings,
             "Heavy atoms": heavy, "QED": round(float(QED.qed(mol)), 3)}
    checks = {
        "Lipinski": (sum((mw > 500, logp > 5, hbd > 5, hba > 10)) <= 1, "MW ≤ 500, cLogP ≤ 5, HBD ≤ 5, HBA ≤ 10 (one violation tolerated)"),
        "Veber": (rotors <= 10 and tpsa <= 140, "Rotatable bonds ≤ 10 and TPSA ≤ 140 Å²"),
        "Ghose": (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 40 <= Descriptors.MolMR(mol) <= 130 and 20 <= mol.GetNumAtoms() <= 70, "MW 160–480, cLogP -0.4–5.6, MR 40–130, atoms 20–70"),
        "Egan": (logp <= 5.88 and tpsa <= 131.6, "cLogP ≤ 5.88 and TPSA ≤ 131.6 Å²"),
        "Muegge": (200 <= mw <= 600 and -2 <= logp <= 5 and tpsa <= 150 and hbd <= 5 and hba <= 10 and rotors <= 15 and rings <= 7, "Common Muegge drug-likeness bounds"),
    }
    alerts = []
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)
    for match in catalog.GetMatches(mol):
        alerts.append(match.GetDescription())
    return {"mol": mol, "properties": props, "rules": checks, "alerts": alerts, "canonical": Chem.MolToSmiles(mol)}

with st.form("analysis"):
    choice = st.selectbox("Try an example molecule", ["Custom input", *EXAMPLES.keys()])
    smiles = st.text_area("SMILES", value=EXAMPLES.get(choice, ""), placeholder="Paste a SMILES string", height=100)
    submitted = st.form_submit_button("Analyze molecule", type="primary", use_container_width=True)

if submitted:
    result = analyze(smiles)
    if "error" in result:
        st.error(result["error"])
    else:
        mol = result["mol"]
        png = Draw.MolToImage(mol, size=(550, 350))
        image_buffer = BytesIO()
        png.save(image_buffer, format="PNG")
        st.subheader("2D structure")
        st.image(image_buffer.getvalue(), use_container_width=True)
        st.caption(f"Canonical SMILES: `{result['canonical']}`")

        st.subheader("Molecular properties")
        prop_items = list(result["properties"].items())
        for start in range(0, len(prop_items), 2):
            cols = st.columns(2)
            for col, (name, value) in zip(cols, prop_items[start:start + 2]):
                col.metric(name, value)

        st.subheader("Drug-likeness rules")
        for name, (passed, explanation) in result["rules"].items():
            st.markdown(f"{'✅ PASS' if passed else '❌ FAIL'} **{name}** — {explanation}")

        st.subheader("PAINS / Brenk alerts")
        if result["alerts"]:
            for alert in result["alerts"]:
                st.warning(alert)
        else:
            st.success("No PAINS or Brenk catalog alerts matched.")

        labels = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"]
        values = [min(mw / lim, 1.5) for mw, lim in zip(
            [result["properties"]["Molecular weight (Da)"], result["properties"]["cLogP"], result["properties"]["TPSA (Å²)"], result["properties"]["HBD"], result["properties"]["HBA"], result["properties"]["Rotatable bonds"]],
            [500, 5.6, 140, 5, 10, 10])]
        chart = go.Figure(go.Scatterpolar(r=values, theta=labels, fill="toself", name="Molecule"))
        chart.add_trace(go.Scatterpolar(r=[1] * len(labels), theta=labels, name="Reference limit", line_dash="dash"))
        chart.update_layout(title="Normalized descriptor comparison", polar=dict(radialaxis=dict(range=[0, 1.5])), height=380)
        st.plotly_chart(chart, use_container_width=True)

st.divider()
st.subheader("Batch CSV")
st.caption("Upload a CSV with a `smiles` column. Free tier limit: 20 rows.")
upload = st.file_uploader("Choose CSV", type=["csv"])
if upload is not None:
    try:
        frame = pd.read_csv(upload)
        if "smiles" not in frame.columns:
            st.error("CSV must contain a column named 'smiles'.")
        elif len(frame) > 20:
            st.error("Free tier supports up to 20 molecules per batch.")
        else:
            output = []
            for idx, value in frame["smiles"].items():
                row = analyze("" if pd.isna(value) else str(value))
                if "error" in row:
                    output.append({"row": idx, "smiles": value, "status": "invalid", "error": row["error"]})
                else:
                    output.append({"row": idx, "smiles": value, "status": "ok", **row["properties"], "error": ""})
            table = pd.DataFrame(output)
            st.dataframe(table, use_container_width=True, hide_index=True)
            st.download_button("Download results CSV", table.to_csv(index=False).encode(), "druglikeness_results.csv", "text/csv", use_container_width=True)
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
        st.error(f"Could not read the CSV: {exc}")

with st.expander("Pro (feature-flag demo)"):
    st.write("Payment integration is TODO. Configure a key in Streamlit secrets to activate this stub.")
    st.text_input("License key", type="password", key="pro_key")
    try:
        expected = st.secrets.get("DRUGLIKENESS_PRO_KEY", "")
    except Exception:
        expected = ""
    if expected and st.session_state.get("pro_key") == expected:
        st.success("Pro key accepted. (PDF and similarity features are not enabled in this demo build.)")
    elif st.session_state.get("pro_key"):
        st.warning("Key not accepted or Pro is not configured.")
