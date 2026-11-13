from __future__ import annotations

from typing import Any

import plotly.graph_objects as go


IDEAL_RANGES: dict[str, tuple[float, float]] = {
    "Molecular weight": (200.0, 500.0),
    "cLogP": (-0.4, 5.6),
    "TPSA": (20.0, 140.0),
    "HBD": (0.0, 5.0),
    "HBA": (0.0, 10.0),
    "Rotatable bonds": (0.0, 10.0),
}


def build_radar_chart(metrics: dict[str, Any]) -> go.Figure:
    """Build a radar chart comparing the molecule to typical drug-like ranges."""
    categories = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"]
    current_values = [
        float(metrics.get("Molecular weight", 0.0)),
        float(metrics.get("cLogP", 0.0)),
        float(metrics.get("TPSA", 0.0)),
        float(metrics.get("HBD", 0.0)),
        float(metrics.get("HBA", 0.0)),
        float(metrics.get("Rotatable bonds", 0.0)),
    ]
    ideal_values = [500.0, 5.6, 140.0, 5.0, 10.0, 10.0]

    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=current_values, theta=categories, fill="toself", name="Molecule", line_color="#2563eb", fillcolor="rgba(37, 99, 235, 0.25)"))
    fig.add_trace(go.Scatterpolar(r=ideal_values, theta=categories, fill="toself", name="Ideal range", line_color="#10b981", fillcolor="rgba(16, 185, 129, 0.1)"))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, max(ideal_values) * 1.1])),
        template="plotly_white",
        height=400,
        margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


def build_bar_chart(metrics: dict[str, Any]) -> go.Figure:
    """Build a bar chart compared against the ideal ranges."""
    properties = [
        ("Molecular weight", "MW", 500.0),
        ("cLogP", "cLogP", 5.6),
        ("TPSA", "TPSA", 140.0),
        ("HBD", "HBD", 5.0),
        ("HBA", "HBA", 10.0),
        ("Rotatable bonds", "RotB", 10.0),
    ]
    labels: list[str] = []
    values: list[float] = []
    limits: list[float] = []
    for name, short_name, ideal_limit in properties:
        labels.append(short_name)
        values.append(float(metrics.get(name, 0.0)))
        limits.append(float(ideal_limit))

    fig = go.Figure(data=[
        go.Bar(name="Molecule", x=labels, y=values, marker_color="#2563eb"),
        go.Bar(name="Ideal limit", x=labels, y=limits, marker_color="#10b981"),
    ])
    fig.update_layout(template="plotly_white", barmode="group", height=400, margin=dict(l=20, r=20, t=20, b=20))
    return fig
