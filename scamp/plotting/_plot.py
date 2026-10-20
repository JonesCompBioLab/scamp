"""
Basic utilities to plot results from scAmp.
"""

import pandas as pd
import plotly.graph_objects as go


def plot_scamp_predictions_plotly(
    prediction_df: pd.DataFrame,
    output_file_name: str,
    title: str = None,
):
    """Plots and saves a Plotly graph of scAmp predictions.

    Saves the file as an interactive HTML file.

    Args:
        prediction_df: A prediction dataframe as returned from a scAmp pipeline.
            Has at least the following columns: `mean`, `var`, `prediction`, and
            `gene`.
        output_file_name: Where to write the output file.
    """
    fig = go.Figure()

    false_predictions = prediction_df[prediction_df["prediction"] == False]
    true_predictions = prediction_df[prediction_df["prediction"] == True]

    trace1 = go.Scatter(
        x=false_predictions["mean"],
        y=false_predictions["var"],
        mode="markers",
        name="non-ecDNA",
        hovertext=list(false_predictions["gene"]),
    )

    trace2 = go.Scatter(
        x=true_predictions["mean"],
        y=true_predictions["var"],
        mode="markers",
        name="ecDNA",
        hovertext=list(true_predictions["gene"]),
    )

    fig.add_trace(trace1)
    fig.add_trace(trace2)

    if title is None:
        title = "scAmp Predictions"

    fig.update_layout(
        title=dict(text=title, font=dict(size=25)),
        xaxis=dict(title=dict(text="Copy-number Mean", font=dict(size=18))),
        yaxis=dict(title=dict(text="Copy-number Variance", font=dict(size=18))),
        legend=dict(font = dict(size = 18, color = "black"))
    )
    fig.write_html(output_file_name)
