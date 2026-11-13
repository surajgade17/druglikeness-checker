# DrugLikeness Checker

A lightweight Streamlit app for screening small molecules by SMILES and comparing them against common oral drug-likeness rules.

## Features

- Analyze a single molecule using SMILES
- View 2D structure, descriptors, rule checks, and alerts
- Compare values against ideal ranges with Plotly charts
- Upload a CSV with a `smiles` column for batch screening
- Use built-in example molecules or paste your own

## Local setup

```bash
python -m venv .venv
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
# macOS/Linux
# source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Notes

This project is intended for research and educational use only, not as a medical or safety prediction tool.
