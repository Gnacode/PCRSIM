import os
import sys
import primer3
from dash import dcc, html, Input, Output, State, dash_table
import dash
from Bio import Entrez, SeqIO
import plotly.graph_objects as go
import matplotlib.colors as mcolors  # For converting named colors to RGB

# Set working directory
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

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

# Function to visualize DNA sequence (horizontal dot plot)
# Function to visualize DNA sequence (horizontal dot plot with dynamic tick adjustment)
def plot_sequence_horizontal_dots(sequence, bases_per_row=50):
    if not sequence:
        fig = go.Figure()
        fig.update_layout(
            title="[no sequence listed]",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    colors = {"A": "red", "T": "green", "C": "blue", "G": "yellow"}
    x_vals, y_vals, dot_colors, hover_texts = [], [], [], []

    for i, base in enumerate(sequence):
        x_vals.append(i // bases_per_row)  # Row index
        y_vals.append(i % bases_per_row)  # Column index within the row
        dot_colors.append(colors.get(base, "grey"))  # Color based on nucleotide
        hover_texts.append(f"Base: {base}<br>Position: {i + 1}")

    trace = go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="markers",
        marker=dict(size=5, color=dot_colors),
        text=hover_texts,
        hoverinfo="text",
    )

    # Calculate the total number of rows for tick adjustment
    total_rows = max(x_vals) + 1
    tick_interval = max(1, total_rows // 10)  # Show 1/10 of x-axis labels when fully zoomed out

    layout = go.Layout(
        title="DNA Sequence Visualization (Rotated Dot Plot)",
        xaxis=dict(
            title="Row Index",
            tickvals=list(range(0, total_rows, tick_interval)),  # Show labels at intervals
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            title="Column Index (Base Position in Row)",
            tickvals=list(range(0, bases_per_row)),  # Show labels for each column
            autorange="reversed",
            showgrid=True,
        ),
        height=800,  # Adjust height dynamically based on rows
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(color="black"),
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

    try:
        results = primer3.bindings.designPrimers(seq_args, global_args)
        primer_pairs = []

        if 'PRIMER_LEFT_NUM_RETURNED' in results and results['PRIMER_LEFT_NUM_RETURNED'] > 0:
            for i in range(results['PRIMER_LEFT_NUM_RETURNED']):
                primer_pairs.append({
                    "Index": i + 1,
                    "Left Primer": results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                    "Right Primer": results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                    "Penalty": results[f'PRIMER_PAIR_{i}_PENALTY'],
                    "Left Position": results[f'PRIMER_LEFT_{i}'][0],
                    "Right Position": results[f'PRIMER_RIGHT_{i}'][0],
                    "Product Size": results[f'PRIMER_RIGHT_{i}'][0] + results[f'PRIMER_RIGHT_{i}'][1] - results[f'PRIMER_LEFT_{i}'][0],
                })

        return primer_pairs

    except Exception as e:
        print(f"Primer design error: {e}")
        return []

# Initialize Dash app
app = dash.Dash(__name__)

# Layout for Dash app
app.layout = html.Div([
    html.H1("PCR Primer Design Dashboard"),
    dcc.Input(id="fasta-id", placeholder="Enter FASTA ID...", type="text", style={"width": "50%"}),
    html.Button("Fetch and Design Primers", id="submit-button", n_clicks=0),
    dcc.Graph(id="sequence-plot", style={"width": "100%"}),  # Full width for the plot
    dash_table.DataTable(
        id="primer-table",
        columns=[{"name": i, "id": i} for i in ["Index", "Left Primer", "Right Primer", "Penalty", "Left Position", "Right Position", "Product Size"]],
        style_table={"overflowX": "auto", "width": "100%"},  # Align with plot width
        style_as_list_view=True,
    ),
])

# Function to visualize DNA sequence with primers
# Function to visualize DNA sequence with primers and equal spacing
import itertools

def plot_sequence_with_primers(sequence, primer_pairs, target_aspect_ratio=1.5):
    if not sequence:
        fig = go.Figure()
        fig.update_layout(
            title="[no sequence listed]",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    # Define colors for primers
    color_palette = itertools.cycle([
        "blue", "red", "green", "orange", "purple", "brown", "cyan", "magenta"
    ])
    primer_colors = [next(color_palette) for _ in range(len(primer_pairs))]

    # Determine optimal bases per row dynamically
    total_bases = len(sequence)
    bases_per_row = int((total_bases * target_aspect_ratio) ** 0.5)
    total_rows = (total_bases + bases_per_row - 1) // bases_per_row

    # Generate sequence data
    x_vals, y_vals, text_vals, hover_texts = [], [], [], []
    for i, base in enumerate(sequence):
        x_vals.append(i % bases_per_row)  # Column index
        y_vals.append(i // bases_per_row)  # Row index
        text_vals.append(base)
        hover_texts.append(f"Base: {base}<br>Position: {i + 1}")

    # Add sequence as text
    trace = go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="text",
        text=text_vals,
        textfont=dict(size=10),
        hoverinfo="text",
        hovertext=hover_texts,
    )

    # Add primers and product lines
    shapes = []
    primer_hover_text = []
    for index, primer in enumerate(primer_pairs):
        primer_color_name = primer_colors[index]
        primer_color_rgb = mcolors.to_rgb(primer_color_name)  # Convert name to RGB
        primer_color_rgba = f"rgba({int(primer_color_rgb[0] * 255)},{int(primer_color_rgb[1] * 255)},{int(primer_color_rgb[2] * 255)},0.2)"

        left_start = primer["Left Position"]
        left_end = left_start + len(primer["Left Primer"]) - 1
        right_start = primer["Right Position"]
        right_end = right_start + len(primer["Right Primer"]) - 1

        # Map primer positions to rows and columns
        left_x0, left_y0 = left_start % bases_per_row, left_start // bases_per_row
        left_x1, left_y1 = left_end % bases_per_row, left_end // bases_per_row
        right_x0, right_y0 = right_start % bases_per_row, right_start // bases_per_row
        right_x1, right_y1 = right_end % bases_per_row, right_end // bases_per_row

        # Add shapes for primers
        shapes.extend([
            # Left primer
            dict(
                type="rect",
                x0=left_x0 - 0.5,
                x1=left_x1 + 0.5,
                y0=left_y0 - 0.5,
                y1=left_y1 + 0.5,
                line=dict(color=primer_color_name, width=2),
                fillcolor=primer_color_rgba,
            ),
            # Right primer
            dict(
                type="rect",
                x0=right_x0 - 0.5,
                x1=right_x1 + 0.5,
                y0=right_y0 - 0.5,
                y1=right_y1 + 0.5,
                line=dict(color=primer_color_name, width=2),
                fillcolor=primer_color_rgba,
            ),
        ])

        # Add a line connecting the primers (product line)
        shapes.append(
            dict(
                type="line",
                x0=left_x1 + 0.5,
                y0=left_y1,
                x1=right_x0 - 0.5,
                y1=right_y0,
                line=dict(color=primer_color_name, width=2, dash="dot"),
            )
        )

        # Primer hover text
        primer_hover_text.append(f"Primer Pair {index + 1}<br>"
                                 f"Left Primer: {primer['Left Primer']}<br>"
                                 f"Right Primer: {primer['Right Primer']}<br>"
                                 f"Product Size: {primer['Product Size']}")

    # Add hover markers for primers
    primer_trace = go.Scatter(
        x=[primer["Left Position"] % bases_per_row for primer in primer_pairs],
        y=[primer["Left Position"] // bases_per_row for primer in primer_pairs],
        mode="markers",
        marker=dict(size=10, color=primer_colors),
        text=primer_hover_text,
        hoverinfo="text",
        name="Primers",
    )

    # Dynamic height and width calculation
    plot_height = 800
    plot_width = int(plot_height * target_aspect_ratio)

    # Define layout
    layout = go.Layout(
        title="DNA Sequence Visualization with Primer Mapping",
        xaxis=dict(
            title="Column Index (Base Position in Row)",
            tickvals=list(range(0, bases_per_row, max(1, bases_per_row // 10))),
            showgrid=True,
            zeroline=False,
            scaleanchor="y",
        ),
        yaxis=dict(
            title="Row Index",
            tickvals=list(range(0, total_rows, max(1, total_rows // 10))),
            autorange="reversed",
            showgrid=True,
        ),
        shapes=shapes,
        height=plot_height,
        width=plot_width,
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(color="black"),
    )

    return go.Figure(data=[trace, primer_trace], layout=layout)



# Update the Dash callback to use this new plot function
@app.callback(
    [Output('sequence-plot', 'figure'),
     Output('primer-table', 'data')],
    [Input('submit-button', 'n_clicks')],
    [State('fasta-id', 'value')]
)
def update_dashboard(n_clicks, fasta_id):
    if not fasta_id:
        print("No FASTA ID provided.")
        return plot_sequence_with_primers(None, []), []  # Placeholder and empty table

    try:
        print(f"Fetching sequence for FASTA ID: {fasta_id}")
        sequence = fetch_sequence_by_fasta(fasta_id)
        print("Sequence fetched successfully.")
        primers = design_primers(sequence)
        print(f"Designed {len(primers)} primer pairs.")
        fig = plot_sequence_with_primers(sequence, primers)
        return fig, primers
    except Exception as e:
        print(f"Error: {e}")
        return plot_sequence_with_primers(None, []), [{"Error": str(e)}]

# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
