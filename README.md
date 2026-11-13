# DrugLikeness Checker

A small Streamlit application for exploratory screening of molecular structures from SMILES. It calculates common molecular descriptors, evaluates five drug-likeness filters, checks a small set of structural-alert SMARTS patterns, and supports CSV batch analysis.

## Features

- Single-molecule analysis with 2D structure image, canonical SMILES, descriptors, filter results, alert matches, and Plotly comparison charts.
- Built-in aspirin, gefitinib, and erlotinib examples.
- Batch CSV analysis with a `smiles` column, a 20-molecule limit, and downloadable results.

## Run locally (Windows PowerShell)

Use Python 3.11 or 3.12. RDKit wheel availability depends on the Python version and platform.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open the local URL printed by Streamlit, usually http://localhost:8501. On macOS/Linux, create the environment with `python3 -m venv .venv` and activate it with `source .venv/bin/activate`.

## Run tests

```bash
python -m pytest -q
```

## Limitations

For research and education only, not medical advice. Rule filters and SMARTS alerts are heuristics; they do not establish a molecule's safety, toxicity, efficacy, or clinical suitability. Some descriptors, including SA score, may be unavailable depending on the RDKit installation.

## Third-party software

The project code is MIT-licensed. RDKit is distributed under a BSD-style license; pandas and Streamlit use BSD-style licenses; Plotly uses the MIT license. Consult upstream packages for full license terms.
