from __future__ import annotations

from typing import Any

import plotly.graph_objects as go

IDEAL_RANGES: dict[str, tuple[float, float]] = {
    "Molecular weight": (200.0, 500.0), "cLogP": (-0.4, 5.6), "TPSA": (20.0, 140.0),
    "HBD": (0.0, 5.0), "HBA": (0.0, 10.0), "Rotatable bonds": (0.0, 10.0),
}


def build_radar_chart(metrics: dict[str, Any]) -> go.Figure:
    """Build radar chart of molecular descriptors and upper-range reference."""
    categories = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"]
    values = [float(metrics.get(k, 0.0) or 0.0) for k in ("Molecular weight", "cLogP", "TPSA", "HBD", "HBA", "Rotatable bonds")]
    limits = [500.0, 5.6, 140.0, 5.0, 10.0, 10.0]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=values, theta=categories, fill="toself", name="Molecule", line_color="#2563eb"))
    fig.add_trace(go.Scatterpolar(r=limits, theta=categories, fill="toself", name="Reference limits", line_color="#10b981", opacity=0.45))
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, max(limits) * 1.1])), template="plotly_white", height=400, margin=dict(l=20, r=20, t=25, b=20))
    return fig


def build_bar_chart(metrics: dict[str, Any]) -> go.Figure:
    """Build grouped descriptor and reference-limit bar chart."""
    properties = [("Molecular weight", "MW", 500.0), ("cLogP", "cLogP", 5.6), ("TPSA", "TPSA", 140.0), ("HBD", "HBD", 5.0), ("HBA", "HBA", 10.0), ("Rotatable bonds", "RotB", 10.0)]
    labels = [short for _, short, _ in properties]
    values = [float(metrics.get(key, 0.0) or 0.0) for key, _, _ in properties]
    limits = [limit for _, _, limit in properties]
    fig = go.Figure(data=[go.Bar(name="Molecule", x=labels, y=values, marker_color="#2563eb"), go.Bar(name="Reference limit", x=labels, y=limits, marker_color="#10b981")])
    fig.update_layout(template="plotly_white", barmode="group", height=400, margin=dict(l=20, r=20, t=25, b=20))
    return fig
