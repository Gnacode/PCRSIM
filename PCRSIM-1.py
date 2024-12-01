import os
import sys
import primer3
from dash import dcc, html, Input, Output, State, dash_table
import dash
from Bio import Entrez, SeqIO
import plotly.graph_objects as go

# Set working directory
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

# Global variables
primer_pairs_simple = []
selected_dna_sequence = ""

# Helper function to clean DNA sequences
def clean_dna_sequence(dna_sequence):
    clean_sequence = dna_sequence.replace(" ", "").upper()
    if not set(clean_sequence).issubset({"A", "T", "C", "G"}):
        raise ValueError("Invalid DNA sequence: must only contain A, T, C, G.")
    return clean_sequence

# Function to fetch DNA sequence using FASTA ID from NCBI
def fetch_sequence_by_fasta(fasta_id):
    Entrez.email = "your_email@example.com"  # Replace with your email
    handle = Entrez.efetch(db="nucleotide", id=fasta_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)

# Function to visualize DNA sequence
def plot_sequence(sequence):
    bases_per_row = 50
    x_vals, y_vals, colors, texts = [], [], [], []
    base_colors = {"A": "lightpink", "T": "lightgreen", "C": "cyan", "G": "lightyellow"}

    for i, base in enumerate(sequence):
        x_vals.append(i % bases_per_row)
        y_vals.append(i // bases_per_row)
        colors.append(base_colors.get(base, "grey"))
        texts.append(f"Base: {base}<br>Position: {i + 1}")

    trace = go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="markers",
        marker=dict(size=5, color=colors),
        text=texts,
        hoverinfo="text",
    )

    layout = go.Layout(
        title="DNA Sequence Visualization",
        xaxis=dict(title="Position in Row", showgrid=False),
        yaxis=dict(title="Row", autorange="reversed"),
    )

    return go.Figure(data=[trace], layout=layout)

# Function to design primers using Primer3
def design_primers(sequence):
    seq_args = {
        'SEQUENCE_ID': 'example',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],
    }

    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_PRODUCT_SIZE_RANGE': [[100, 300]],
        'PRIMER_NUM_RETURN': 10,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }

    results = primer3.designPrimers(seq_args, global_args)
    primer_pairs = []

    if 'PRIMER_LEFT_NUM_RETURNED' in results and results['PRIMER_LEFT_NUM_RETURNED'] > 0:
        for i in range(results['PRIMER_LEFT_NUM_RETURNED']):
            primer_pairs.append({
                "Left Primer": results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                "Right Primer": results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                "Penalty": results[f'PRIMER_PAIR_{i}_PENALTY'],
                "Left Position": results[f'PRIMER_LEFT_{i}'][0],
                "Right Position": results[f'PRIMER_RIGHT_{i}'][0],
                "Product Size": results[f'PRIMER_RIGHT_{i}'][0] + results[f'PRIMER_RIGHT_{i}'][1] - results[f'PRIMER_LEFT_{i}'][0],
            })

    return primer_pairs

# Initialize Dash app
app = dash.Dash(__name__)

# Layout for Dash app
app.layout = html.Div([
    html.H1("PCR Primer Design Dashboard"),
    dcc.Input(id="fasta-id", placeholder="Enter FASTA ID...", type="text", style={"width": "50%"}),
    html.Button("Fetch and Design Primers", id="submit-button", n_clicks=0),
    dcc.Graph(id="sequence-plot"),
    dash_table.DataTable(
        id="primer-table",
        columns=[
            {"name": "Left Primer", "id": "Left Primer"},
            {"name": "Right Primer", "id": "Right Primer"},
            {"name": "Penalty", "id": "Penalty"},
            {"name": "Left Position", "id": "Left Position"},
            {"name": "Right Position", "id": "Right Position"},
            {"name": "Product Size", "id": "Product Size"},
        ],
        style_table={"overflowX": "auto"},
    ),
])

# Callbacks for interactivity
@app.callback(
    [Output("sequence-plot", "figure"),
     Output("primer-table", "data")],
    [Input("submit-button", "n_clicks")],
    [State("fasta-id", "value")],
)
def update_dashboard(n_clicks, fasta_id):
    if n_clicks == 0 or not fasta_id:
        return {}, []

    try:
        sequence = fetch_sequence_by_fasta(fasta_id)
        fig = plot_sequence(sequence)
        primers = design_primers(sequence)
        return fig, primers
    except Exception as e:
        return {}, [{"Error": str(e)}]

# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
