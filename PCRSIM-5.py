import os
import sys
import primer3
from dash import dcc, html, Input, Output, State, dash_table
import dash
from Bio import Entrez, SeqIO
import plotly.graph_objects as go
import matplotlib.colors as mcolors  # For converting named colors to RGB
import itertools  # For cycling through colors

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
                    "Product Size": results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                })

        return primer_pairs

    except Exception as e:
        print(f"Primer design error: {e}")
        return []

# Function to add product region shapes
def add_product_region_shapes(shapes, product_start, product_end, bases_per_row, light_color):
    start_row = product_start // bases_per_row
    end_row = product_end // bases_per_row
    for row in range(start_row, end_row + 1):
        row_start_pos = row * bases_per_row
        row_end_pos = row_start_pos + bases_per_row - 1

        current_start = max(product_start, row_start_pos)
        current_end = min(product_end, row_end_pos)

        x0 = current_start % bases_per_row
        y0 = row

        x1 = current_end % bases_per_row
        y1 = row

        shapes.append(
            dict(
                type="rect",
                x0=x0 - 0.5,
                x1=x1 + 0.5,
                y0=y0 - 0.5,
                y1=y1 + 0.5,
                line=dict(color=light_color, width=1),
                fillcolor=light_color,
            )
        )

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

    # Define colors for bases and primers
    base_colors = {"A": "red", "T": "green", "C": "blue", "G": "yellow"}
    color_palette = itertools.cycle([
        "blue", "red", "green", "orange", "purple", "brown", "cyan", "magenta"
    ])
    primer_colors = [next(color_palette) for _ in range(len(primer_pairs))]

    # Determine optimal bases per row dynamically
    total_bases = len(sequence)
    bases_per_row = int((total_bases * target_aspect_ratio) ** 0.5)
    total_rows = (total_bases + bases_per_row - 1) // bases_per_row

    # Generate sequence data
    x_vals, y_vals, dot_colors, text_vals, hover_texts = [], [], [], [], []
    for i, base in enumerate(sequence):
        x_vals.append(i % bases_per_row)  # Column index
        y_vals.append(i // bases_per_row)  # Row index
        dot_colors.append(base_colors.get(base, "gray"))  # Use base colors
        text_vals.append(base)  # Base letters for overlay
        hover_texts.append(f"Base: {base}<br>Position: {i + 1}")

    # Create traces for dots and letters
    dots_trace = go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="markers",
        marker=dict(
            size=6,
            color=[f"rgba({int(mcolors.to_rgb(base_colors[base])[0]*255)},"
                   f"{int(mcolors.to_rgb(base_colors[base])[1]*255)},"
                   f"{int(mcolors.to_rgb(base_colors[base])[2]*255)},0.2)" for base in sequence],
        ),
        name="Colored Dots",  # Legend for dots
        hoverinfo="text",
        hovertext=hover_texts,
        visible=True,  # Default ON
    )

    letters_trace = go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="text",
        text=text_vals,
        textfont=dict(size=5),
        hoverinfo="text",
        hovertext=hover_texts,
        name="Base Letters",  # Legend for letters
        visible=False,  # Default OFF
    )

    # Add individual primer and product traces
    primers_traces = []
    products_traces = {}
    for index, primer in enumerate(primer_pairs):
        primer_color_name = primer_colors[index]
        primer_color_rgb = mcolors.to_rgb(primer_color_name)  # Convert name to RGB
        dark_color = f"rgba({int(primer_color_rgb[0]*255)},{int(primer_color_rgb[1]*255)},{int(primer_color_rgb[2]*255)},0.8)"
        light_color = f"rgba({int(primer_color_rgb[0]*255)},{int(primer_color_rgb[1]*255)},{int(primer_color_rgb[2]*255)},0.3)"

        left_start = primer["Left Position"]
        left_end = left_start + len(primer["Left Primer"]) - 1
        right_start = primer["Right Position"]
        right_end = right_start + len(primer["Right Primer"]) - 1

        # Map primer positions to rows and columns
        left_x0, left_y0 = left_start % bases_per_row, left_start // bases_per_row
        left_x1, left_y1 = left_end % bases_per_row, left_end // bases_per_row
        right_x0, right_y0 = right_start % bases_per_row, right_start // bases_per_row
        right_x1, right_y1 = right_end % bases_per_row, right_end // bases_per_row

        # Primer rectangles
        primers_traces.append(go.Scatter(
            x=[left_x0, left_x1, None, right_x0, right_x1],
            y=[left_y0, left_y1, None, right_y0, right_y1],
            mode="lines",
            line=dict(color=dark_color, width=6),
            name=f"Primer Pair {index + 1}",
            hoverinfo="text",
        ))

        # Product regions
        product_start = left_end + 1
        product_end = right_start - 1

        if product_start <= product_end:
            product_trace = products_traces.get(index + 1, [])
            start_row = product_start // bases_per_row
            end_row = product_end // bases_per_row
            for row in range(start_row, end_row + 1):
                row_start_pos = row * bases_per_row
                row_end_pos = row_start_pos + bases_per_row - 1

                current_start = max(product_start, row_start_pos)
                current_end = min(product_end, row_end_pos)

                x0 = current_start % bases_per_row
                y0 = row

                x1 = current_end % bases_per_row
                y1 = row

                product_trace.append((x0, y0, x1, y1))

            products_traces[index + 1] = product_trace

    # Consolidate product traces
    product_scatter_traces = []
    for product_id, fragments in products_traces.items():
        x_vals, y_vals = [], []
        for x0, y0, x1, y1 in fragments:
            x_vals += [x0, x1, None]
            y_vals += [y0, y1, None]

        product_scatter_traces.append(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="lines",
            line=dict(color=light_color, width=3),
            name=f"Product {product_id}",
        ))

    # Combine traces
    traces = [dots_trace, letters_trace] + primers_traces + product_scatter_traces

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
        height=plot_height,
        width=plot_width,
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(color="black"),
        legend=dict(
            itemsizing="constant",
            title="Legend",
        ),
    )

    return go.Figure(data=traces, layout=layout)










# Initialize Dash app
app = dash.Dash(__name__)

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

@app.callback(
    [Output('sequence-plot', 'figure'),
     Output('primer-table', 'data')],
    [Input('submit-button', 'n_clicks')],
    [State('fasta-id', 'value')]
)
def update_dashboard(n_clicks, fasta_id):
    if not fasta_id:
        return plot_sequence_with_primers(None, []), []

    try:
        sequence = fetch_sequence_by_fasta(fasta_id)
        primers = design_primers(sequence)
        fig = plot_sequence_with_primers(sequence, primers)
        return fig, primers
    except Exception as e:
        return plot_sequence_with_primers(None, []), [{"Error": str(e)}]

if __name__ == "__main__":
    app.run_server(debug=True)
