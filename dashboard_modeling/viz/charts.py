"""Plotly charts for the dashboard."""
import plotly.express as px


def strength_histogram(edges, title="Edge strength distribution"):
    """Histogram of bootstrap edge strengths (all scored edges)."""
    fig = px.histogram(edges, x="strength", nbins=40)
    fig.update_layout(
        height=320, margin=dict(t=40, b=10),
        title=title, xaxis_title="Edge strength", yaxis_title="Count",
    )
    fig.add_vline(x=0.85, line_dash="dash", line_color="#C44E52")
    return fig
