"""Plotly charts for the dashboard."""
import numpy as np
import plotly.graph_objects as go

# Default edge-strength threshold, shared by the Explorer slider and the
# histogram reference line so the two never drift apart.
DEFAULT_STRENGTH_THRESHOLD = 0.85


def strength_histogram(edges, title="Edge strength distribution"):
    """Histogram of bootstrap edge strengths (all scored edges).

    Pre-aggregates server-side (40 bars) rather than shipping every raw
    strength value to the browser — matters for L_all (~123k edges), whose
    figure would otherwise be rebuilt and re-transmitted on every rerun.
    """
    counts, bin_edges = np.histogram(
        edges["strength"].dropna(), bins=40, range=(0.0, 1.0)
    )
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    fig = go.Figure(go.Bar(x=centers, y=counts, width=1 / 40,
                           marker_color="#4C72B0"))
    fig.update_layout(
        height=320, margin=dict(t=40, b=10), bargap=0,
        title=title, xaxis_title="Edge strength", yaxis_title="Count",
    )
    fig.add_vline(x=DEFAULT_STRENGTH_THRESHOLD, line_dash="dash",
                  line_color="#C44E52")
    return fig
