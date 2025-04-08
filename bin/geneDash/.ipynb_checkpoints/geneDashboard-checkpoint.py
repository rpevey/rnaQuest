import dash
from dash import dcc, html, Input, Output, State
import dash.dash_table
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

# Load the merged data
df = pd.read_csv("dash_expression_data.csv", na_values=["NA", "NaN", ""])

# Clean for safety
df['gene'] = df['gene'].astype(str)
df['condition'] = df['condition'].astype(str)
# Convert to numeric, ensuring float64 dtype
df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce", downcast=None)
df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce', downcast=None)
#df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")

# Set up gene list
genes = sorted(df['gene'].unique())

# Dash app setup
app = dash.Dash(__name__)

default_gene = df.loc[df['log2FoldChange'].abs().idxmax(), 'gene']

#Viridis color palette
viridis_colors = sample_colorscale("Viridis", [0.75, 0.25])

app.layout = html.Div([
  dcc.Dropdown(
    id="gene-dropdown",
    options=[{"label": g, "value": g} for g in sorted(df["gene"].dropna().unique())],
    # Set a specific gene of interest as the default gene
    #value="LGR6",
    value=default_gene,
    style={"width": "200px"}
  ),
  dcc.Graph(id="expression-plot",
           style={"width": "100%",
                  "margin": "0 auto"}
           ),
  dcc.Checklist(
    id='sig-only-toggle',
    options=[{'label': ' Show only significant genes (FDR < 0.05)', 'value': 'sig'}],
    value=[],  # empty = show all by default
    style={"margin": "10px"}
),
  html.Button("Download graph data (CSV)", id="download-btn", n_clicks=0, style={"margin": "20px"}),
  dcc.Download(id="download-data"),
  html.Div(id="de-info-text", style={"margin": "10px", "textAlign": "center", "fontSize": "14px"}),
  dash.dash_table.DataTable(
    id='expression-table',
    style_table={'overflowX': 'auto', 'margin': '20px auto', 'width': '80%'},
    style_cell={'textAlign': 'left',
                'padding': '5px',
                'fontFamily': 'sans-serif',
                'fontSize': '14px'
               },
    style_data_conditional=[
        {
            'if': {'column_id': 'expression'},
            'textAlign': 'right',
            'fontFamily': 'monospace'
        }
    ],
    style_header={'fontWeight': 'bold'},
    page_size=10
  )
])

@app.callback(
  Output("gene-dropdown", "options"),
  Output("gene-dropdown", "value"),
  Input("sig-only-toggle", "value")
)

def update_gene_list(toggle_values):
    if "sig" in toggle_values:
        filtered = df[df["FDR"] < 0.05]
    else:
        filtered = df

    # Remove NAs and duplicates
    gene_list = sorted(filtered["gene"].dropna().unique())

    # Default gene: max abs log2FC
    if len(gene_list) == 0:
        return [], None
    default_gene = (
        filtered.loc[filtered["log2FoldChange"].abs().idxmax(), "gene"]
        if "log2FoldChange" in filtered.columns else gene_list[0]
    )

    options = [{"label": g, "value": g} for g in gene_list]
    
    return options, default_gene

@app.callback(
  Output("expression-plot", "figure"),
  Input("gene-dropdown", "value")
)

def update_plot(selected_gene):
    #Define dataset for the plot
    data = df[df["gene"] == selected_gene].copy()

     # Log-transform raw counts for plotting
    data["log_counts"] = np.log10(data["counts"] + 1)

    # --- Calculate y-axis limits and ticks ---
    # Ensure data has valid counts
    if data["counts"].notnull().any() and (data["counts"] > 0).any():
        y_max = data["counts"].max()
    else:
        y_max = 10  # fallback if no valid counts are found
    
    max_tick = int(np.ceil(np.log10(y_max + 1)))
    tick_vals = [10**i for i in range(0, max_tick + 1)]

   #Define FDR and LFC text
    padj = data["FDR"].iloc[0]
    lfc = data["log2FoldChange"].iloc[0]
    padj_text = f"FDR = {padj:.2e}" if padj<0.05 else "Not DE"

    fig = go.Figure()

    #Convert conditions from catagorical to numeric to plot boxplots and jittered datapoints on the same field
    conditions = ["TL", "TLE"]  # TL = Control, TLE = Epilepsy
    cond_to_x = {cond: i for i, cond in enumerate(conditions)}
    #Add condition labels.
    cond_labels = {
    "TL": "Control",
    "TLE": "Epilepsy"
    }
    
    for i, cond in enumerate(conditions):
        subset = data[data["condition"] == cond]
        x_val = cond_to_x[cond]
        
        #Boxplot
        fig.add_trace(go.Box(
            y=subset["log_counts"],
            x=[x_val] * len(subset),
            name=cond,
            marker_color="white", #fill color
            fillcolor="white", # explicitly set
            line=dict(color="black"), #line color
            width=0.25,
            boxpoints="outliers",
            showlegend=False,
            hoverinfo="skip"
        ))

        # Dot plot (viridis-colored with jitter)
        x_jittered = [x_val + np.random.uniform(-0.15, 0.15) for _ in range(len(subset))]

        #Scatter points
        fig.add_trace(go.Scatter(
            y=subset["log_counts"],
            x=x_jittered,
            mode="markers",
            marker=dict(
                size=12, 
                color=viridis_colors[i % len(viridis_colors)]
            ),
            hovertext=[
                f"Sample: {s}<br>Counts: {c:.3g}"
                for s, c in zip(subset["sample"], subset["counts"])
            ],
            hoverinfo="text",
            showlegend=False
        ))

        #Add significance bar for DE genes
        # Check if gene is significant
        fdr = data["FDR"].iloc[0]
        lfc = data["log2FoldChange"].iloc[0]
        is_significant = pd.notnull(fdr) and fdr < 0.05

        if is_significant:
            bar_y = np.log10(y_max * 1.5 + 1)
        
            # Horizontal line
            fig.add_trace(go.Scatter(
                x=[0, 1],
                y=[bar_y, bar_y],
                mode="lines",
                line=dict(color="black", width=1),
                showlegend=False,
                hoverinfo="skip"
            ))
        
            # Vertical line 0
            fig.add_trace(go.Scatter(
                x=[0, 0],
                y=[(bar_y*0.98), bar_y],
                mode="lines",
                line=dict(color="black", width=1),
                showlegend=False,
                hoverinfo="skip"
            ))

            # Vertical line 1
            fig.add_trace(go.Scatter(
                x=[1, 1],
                y=[(bar_y*0.98), bar_y],
                mode="lines",
                line=dict(color="black", width=1),
                showlegend=False,
                hoverinfo="skip"
            ))

    fig.update_layout(
        title=f"<i>{selected_gene}</i>",
        yaxis_title="Counts",
        title_x=0.5,
        font=dict(size=12),
        plot_bgcolor="white",
        xaxis=dict(
            range=[-0.5, 1.5],  # tighter than default for 2 categories
        ),
        yaxis=dict(
            type="linear", # Log transformed data but linear axis
            autorange=False,       # Force range to apply
            range=[0, np.ceil(data["log_counts"].max()) * 1.1],
            tickvals=[np.log10(x + 1) for x in tick_vals],  # Transformed for log axis
            ticktext=[str(x) for x in tick_vals],           # Displayed as raw counts
            tickformat=".0f"
        ),
        annotations=[
            dict(
                text=padj_text,
                x=0.5,
                y=np.log10(y_max * 1.9 + 1),
                xref="x",
                yref="y",
                showarrow=False,
                font=dict(size=12)
            )
        ]
    )

    fig.update_yaxes(type="linear", linecolor="black")
    fig.update_xaxes(
        linecolor="black",
        tickangle=-45,
        tickmode="array",
        tickvals=list(cond_to_x.values()),
        ticktext=[cond_labels[cond] for cond in conditions],
        type="linear"
    )

    return fig

@app.callback(
    Output("download-data", "data"),
    Input("download-btn", "n_clicks"),
    State("gene-dropdown", "value"),
    prevent_initial_call=True
)
def download_table(n_clicks, selected_gene):
    if selected_gene is None:
        return dash.no_update

    data_subset = df[df["gene"] == selected_gene].copy()

    columns_to_export = ["sample", "condition", "counts"]
    if "FDR" in data_subset.columns and "log2FoldChange" in data_subset.columns:
        columns_to_export += ["log2FoldChange", "FDR"]

    return dcc.send_data_frame(
        data_subset[columns_to_export].to_csv,
        filename=f"{selected_gene}_expression_table.csv",
        index=False
    )

@app.callback(
    Output("de-info-text", "children"),
    Input("gene-dropdown", "value")
)
def update_de_info(selected_gene):
    if selected_gene is None:
        return ""

    gene_data = df[df["gene"] == selected_gene]

    if "FDR" in gene_data.columns and "log2FoldChange" in gene_data.columns:
        fdr_val = gene_data["FDR"].iloc[0]
        lfc_val = gene_data["log2FoldChange"].iloc[0]

        if pd.notnull(fdr_val) and pd.notnull(lfc_val):
            return f"log₂FC: {lfc_val:.3f} | FDR: {fdr_val:.2e}"
        else:
            return "No differential expression data available for this gene."
    else:
        return "DE information not available."


@app.callback(
    Output("expression-table", "data"),
    Output("expression-table", "columns"),
    Input("gene-dropdown", "value")
)
def update_table(selected_gene):
    if selected_gene is None:
        return [], []

    data_subset = df[df["gene"] == selected_gene].copy()
    #Round expression data to 3 significant digits
    data_subset["counts"] = data_subset["counts"].apply(lambda x: f"{x:.4g}")

    #Choose what to show in the table
    columns_to_show = ["sample", "condition", "counts"]
    #Replace columns for readability
    column_headers = {
    "sample": "Sample",
    "condition": "Condition",
    "counts": "Normalized gene counts"
    }
    
    table_data = data_subset[columns_to_show].to_dict("records")
    table_columns = [{"name": column_headers[col], "id": col} for col in columns_to_show]

    return table_data, table_columns


if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8080)