# DrugLikeness Checker

A simple Streamlit application to explore molecular descriptors and drug-likeness filters using SMILES. All demo features are available without a license key.

## Run locally

Requires Python 3.11 or 3.12 recommended (RDKit wheel support can vary by Python version).

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open the printed local URL (usually http://localhost:8501). On macOS/Linux, create and activate with `python3 -m venv .venv` and `source .venv/bin/activate`.

## Codespaces

Use Code → Codespaces → Create codespace on main. In its terminal run `python -m pip install -r requirements.txt`, then `python -m streamlit run app.py --server.address 0.0.0.0`. Open forwarded port 8501.

## Features and use

- Paste SMILES or select aspirin, gefitinib, or erlotinib, then click Analyze.
- View 2D structure, descriptors, rule checks, alerts, and a Plotly chart.
- Upload a CSV with a `smiles` column; demo batch size is capped at 5,000 rows and results are downloadable.
- No license key is required; payment integration is not implemented.

## Limitations

For research and educational use only. Not medical advice. Descriptors and rule filters are heuristics, not efficacy, toxicity, or clinical predictions. SA score is optional and may show unavailable if RDKit Contrib is absent. Structural alerts need expert interpretation.

## License

Application code is MIT licensed. RDKit is BSD-licensed; Streamlit and pandas are BSD-licensed; Plotly is MIT-licensed. Refer to upstream projects for complete third-party notices.
