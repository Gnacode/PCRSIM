import os
import sys
import primer3
from dash import Input, Output, dcc, html, State, dash_table
import dash_bootstrap_components as dbc
import dash
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import itertools
import numpy as np

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
        results = primer3.bindings.design_primers(seq_args, global_args)
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

def plot_sequence_with_primers(sequence, primer_pairs, target_aspect_ratio=1.5, max_sequence_length=50000):
    """
    Plot a DNA sequence with primers. If the sequence length exceeds `max_sequence_length`,
    only plot primers and products without the full sequence.
    """
    if len(sequence) > max_sequence_length:
        # Skip sequence visualization, show primers and products only
        fig = go.Figure()
        fig.update_layout(
            title=f"DNA Sequence Visualization (sequence length: {len(sequence)} bases, showing primers and products only)",
            xaxis=dict(title="Base Position", showgrid=True),
            yaxis=dict(title="Layer", visible=False),  # Hide unused y-axis
            template="plotly_white"
        )

        for index, primer in enumerate(primer_pairs):
            primer_color = f"rgba({(index * 35) % 255}, {(index * 70) % 255}, {(index * 105) % 255}, 0.8)"
            product_color = f"rgba({(index * 35) % 255}, {(index * 70) % 255}, {(index * 105) % 255}, 0.3)"

            # Primer positions
            left_start = primer["Left Position"]
            left_end = left_start + len(primer["Left Primer"]) - 1
            right_start = primer["Right Position"]
            right_end = right_start + len(primer["Right Primer"]) - 1

            fig.add_trace(go.Scatter(
                x=[left_start, left_end],
                y=[0, 0],
                mode="lines",
                line=dict(color=primer_color, width=4),
                name=f"Primer Pair {index + 1} (Left)"
            ))

            fig.add_trace(go.Scatter(
                x=[right_start, right_end],
                y=[0, 0],
                mode="lines",
                line=dict(color=primer_color, width=4),
                name=f"Primer Pair {index + 1} (Right)"
            ))

            # Product regions
            product_start = left_end + 1
            product_end = right_start - 1
            if product_start <= product_end:
                fig.add_trace(go.Scatter(
                    x=[product_start, product_end],
                    y=[0, 0],
                    mode="lines",
                    line=dict(color=product_color, width=2, dash="dash"),
                    name=f"Product {index + 1}"
                ))

        return fig

    else:

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

        # Consolidate primers and products
        primers_traces = []
        products_traces = []

        # Updated portion of the `plot_sequence_with_primers` function to handle multi-row primers

        for index, primer in enumerate(primer_pairs):
            primer_color_name = primer_colors[index]
            primer_color_rgb = mcolors.to_rgb(primer_color_name)  # Convert name to RGB
            dark_color = f"rgba({int(primer_color_rgb[0]*255)},{int(primer_color_rgb[1]*255)},{int(primer_color_rgb[2]*255)},0.8)"
            light_color = f"rgba({int(primer_color_rgb[0]*255)},{int(primer_color_rgb[1]*255)},{int(primer_color_rgb[2]*255)},0.3)"  # Light shade

            # Fully span the left and right primers
            left_start = primer["Left Position"]
            left_end = left_start + len(primer["Left Primer"]) - 1
            right_start = primer["Right Position"]
            right_end = right_start + len(primer["Right Primer"]) - 1

            # Handle multi-row left primer
            left_x_vals, left_y_vals = [], []
            for pos in range(left_start, left_end + 1):
                left_x_vals.append(pos % bases_per_row)
                left_y_vals.append(pos // bases_per_row)
                if pos % bases_per_row == bases_per_row - 1:
                    left_x_vals.append(None)  # Break line for new row
                    left_y_vals.append(None)

            # Handle multi-row right primer
            right_x_vals, right_y_vals = [], []
            for pos in range(right_start, right_end + 1):
                right_x_vals.append(pos % bases_per_row)
                right_y_vals.append(pos // bases_per_row)
                if pos % bases_per_row == bases_per_row - 1:
                    right_x_vals.append(None)  # Break line for new row
                    right_y_vals.append(None)

            # Add left and right primers as one trace
            primers_traces.append(go.Scatter(
                x=left_x_vals + [None] + right_x_vals,
                y=left_y_vals + [None] + right_y_vals,
                mode="lines",
                line=dict(color=dark_color, width=4),
                name=f"Primer Pair {index + 1}",
                hoverinfo="text",
                hovertext=f"Primer Pair {index + 1}<br>Left Primer: {primer['Left Primer']}<br>Right Primer: {primer['Right Primer']}",
            ))

            # Product regions
            product_start = left_end + 1
            product_end = right_start - 1

            if product_start <= product_end:
                product_x_vals = []
                product_y_vals = []

                start_row = product_start // bases_per_row
                end_row = product_end // bases_per_row

                for row in range(start_row, end_row + 1):
                    row_start = row * bases_per_row
                    row_end = row_start + bases_per_row - 1

                    current_start = max(product_start, row_start)
                    current_end = min(product_end, row_end)

                    x_start = current_start % bases_per_row
                    y_start = current_start // bases_per_row
                    x_end = current_end % bases_per_row
                    y_end = current_end // bases_per_row

                    product_x_vals += [x_start, x_end, None]
                    product_y_vals += [y_start, y_end, None]

                products_traces.append(go.Scatter(
                    x=product_x_vals,
                    y=product_y_vals,
                    mode="lines",
                    line=dict(color=light_color, width=2),
                    name=f"Product {index + 1}",
                    hoverinfo="text",
                    hovertext=f"Product {index + 1}<br>Size: {primer['Product Size']} bases",
                ))


        # Combine traces
        traces = [dots_trace, letters_trace] + primers_traces + products_traces

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


## PCR AMPLICATON SECTON

def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence to be reverse complemented.
    
    Returns:
        str: Reverse complement of the input sequence.
    """
    return str(Seq(sequence).reverse_complement())

# Enthalpy and entropy tables for Nearest-Neighbor model
enthalpy = {
    'AA': -7.9, 'TT': -7.9,
    'AT': -7.2, 'TA': -7.2,
    'CA': -8.5, 'TG': -8.5,
    'GT': -8.4, 'AC': -8.4,
    'CT': -7.8, 'AG': -7.8,
    'GA': -8.2, 'TC': -8.2,
    'GG': -8.0, 'CC': -8.0,
    'GC': -9.8, 'CG': -9.8
}

entropy = {
    'AA': -22.2, 'TT': -22.2,
    'AT': -20.4, 'TA': -20.4,
    'CA': -22.7, 'TG': -22.7,
    'GT': -22.4, 'AC': -22.4,
    'CT': -21.0, 'AG': -21.0,
    'GA': -22.2, 'TC': -22.2,
    'GG': -19.9, 'CC': -19.9,
    'GC': -24.4, 'CG': -24.4
}

def sigmoid(x, x0=0.5, k=1.0):
    return 1 / (1 + np.exp(-k * (x - x0)))

def smooth_curve(points, factor=0.9):
    smoothed_points = []
    for point in points:
        if smoothed_points:
            previous = smoothed_points[-1]
            smoothed_points.append(previous * factor + point * (1 - factor))
        else:
            smoothed_points.append(point)
    return smoothed_points

def calculate_tm(sequence):
    delta_H = 0
    delta_S = 0

    for i in range(len(sequence) - 1):
        pair = sequence[i:i + 2]
        if pair in enthalpy:
            delta_H += enthalpy[pair]
            delta_S += entropy[pair]

    R = 1.987  # cal/mol·K
    C = 50e-9  # DNA duplex concentration in M

    tm = (delta_H * 1000) / (delta_S + R * np.log(C / 4)) - 273.15
    return tm


def adjust_efficiency(original_efficiency, product_length):
    efficiency_drop_per_100bp = 0.15
    adjusted_efficiency = original_efficiency * (1 - (product_length // 100) * efficiency_drop_per_100bp)
    return max(adjusted_efficiency, 0)

def calculate_modified_efficiency(primer_pair, product_length, results, index):
    """
    Calculate a modified amplification efficiency based on primer penalties, product length,
    melting temperature, and GC content.
    """
    left_primer_sequence = primer_pair['Left Primer']
    right_primer_sequence = primer_pair['Right Primer']

    # Calculate melting temperatures for both primers
    left_tm = primer3.calc_tm(left_primer_sequence)
    right_tm = primer3.calc_tm(right_primer_sequence)
    average_tm = (left_tm + right_tm) / 2

    # GC content
    def gc_content(seq):
        return (seq.count("G") + seq.count("C")) / len(seq) * 100

    left_gc = gc_content(left_primer_sequence)
    right_gc = gc_content(right_primer_sequence)
    average_gc = (left_gc + right_gc) / 2

    # Retrieve penalties from the results
    left_penalty = results[f'PRIMER_LEFT_{index}_PENALTY']
    right_penalty = results[f'PRIMER_RIGHT_{index}_PENALTY']

    # Base efficiency considering length
    base_efficiency = adjust_efficiency(1.0, product_length)

    # Penalty for Tm deviation (ideal range: 57°C–63°C)
    tm_deviation_penalty = max(0, abs(average_tm - 60) / 10)

    # Penalty for GC content deviation (ideal range: 40%–60%)
    gc_deviation_penalty = max(0, abs(average_gc - 50) / 20)

    # Combine penalties
    penalty_factor = 1 - ((left_penalty + right_penalty) / 2)
    tm_factor = 1 - tm_deviation_penalty
    gc_factor = 1 - gc_deviation_penalty

    # Combine factors into a final efficiency
    modified_efficiency = base_efficiency * penalty_factor * tm_factor * gc_factor
    return max(modified_efficiency, 0)  # Ensure efficiency is non-negative



def find_product_length(primer_pair, sequence):
    forward_primer, reverse_primer = primer_pair
    reverse_primer = reverse_complement(reverse_primer)

    forward_start = sequence.find(forward_primer)
    if forward_start == -1:
        raise ValueError(f"Forward primer {forward_primer} not found in the sequence.")

    reverse_start = sequence.rfind(reverse_primer)
    if reverse_start == -1:
        raise ValueError(f"Reverse primer {reverse_primer} not found in the sequence.")

    product_length = reverse_start + len(reverse_primer) - forward_start
    if product_length <= 0:
        raise ValueError(f"Invalid product length: {product_length}.")

    return product_length




def simulate_pcr_amplification(primer_pairs, sequence, results):
    amplification_curves = []
    melting_curves = []
    melting_temps = []

    initial_concentration = 1e6
    max_primers = 1e12
    num_cycles = 55

    for index, primer_pair in enumerate(primer_pairs):
        left_primer = primer_pair['Left Primer']
        right_primer = primer_pair['Right Primer']

        # Calculate product length
        product_length = find_product_length((left_primer, right_primer), sequence)

        # Calculate modified efficiency
        efficiency_base = calculate_modified_efficiency(primer_pair, product_length, results, index)

        # PCR amplification simulation
        amplification = [initial_concentration]
        for cycle in range(1, num_cycles):
            efficiency = sigmoid((max_primers - amplification[-1]) / max_primers, x0=0.5, k=5)
            amplification.append(amplification[-1] * (1 + efficiency * efficiency_base))
            amplification[-1] = min(amplification[-1], max_primers)

        amplification_curves.append({
            'x': list(range(num_cycles)),
            'y': smooth_curve(amplification),
            'name': f"Primer Pair {index + 1}"
        })

        # Melting curve simulation
        product_sequence = sequence[sequence.find(left_primer):sequence.rfind(reverse_complement(right_primer)) + len(reverse_complement(right_primer))]
        melting_temp = calculate_tm(product_sequence)
        melting_temps.append(melting_temp)

        temp_range = np.linspace(40, 100, 500)
        fluorescence = -np.gradient(sigmoid(temp_range, x0=melting_temp, k=10))
        melting_curves.append({
            'x': temp_range.tolist(),
            'y': fluorescence.tolist(),
            'name': f"Primer Pair {index + 1}"
        })

    return amplification_curves, melting_curves

## CAPILLARY ELECTROPHORESIS SIMULATION

def simulate_capillary_electrophoresis(amplicon_lengths):
    """
    Simulate Capillary Electrophoresis for given amplicon lengths.
    """
    np.random.seed(0)  # For reproducibility

    # Create a function to generate a Gaussian peak
    def gaussian(x, mu, sigma, amplitude):
        return amplitude * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    # Ladder positions (from 100 bp to 1000 bp)
    ladder_bp = np.arange(100, 1100, 100)

    # Create a range of x values for plotting
    x = np.linspace(0, 1200, 10000)

    # Generate ladder and sample peaks
    ladder = np.zeros_like(x)

    ladder_amplitude = 5e5

    for bp in ladder_bp:
        ladder += gaussian(x, bp, 10 - 0.005 * bp, ladder_amplitude)  # Inverted spread with distance
        noise = 0.001 * ladder_amplitude  # 10% noise
        ladder += np.random.normal(0, noise, x.shape)

    # Generate amplicon peaks
    amplicon_peaks = np.zeros_like(x)
    sample_amplitude = 1e6  # Adjust as needed
    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]

    # Create Plotly figure
    fig = go.Figure()

    # Add ladder data
    fig.add_trace(
        go.Scatter(
            x=x,
            y=ladder,
            mode="lines",
            name="Ladder",
            line=dict(color="red", width=1),
        )
    )

    # Add amplicon data
    for i, bp in enumerate(amplicon_lengths):
        amplicon_peak = gaussian(x, bp, 10 - 0.005 * bp, sample_amplitude)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=amplicon_peak,
                mode="lines",
                name=f"Amplicon {bp}bp",
                line=dict(color=colors[i % len(colors)], width=1),
            )
        )

    # Update layout properties
    fig.update_layout(
        title="Simulated Capillary Electrophoresis",
        xaxis_title="Migration Distance (bp)",
        yaxis_title="Amplicon Concentration",
        xaxis=dict(autorange="reversed"),  # Invert x-axis
        template="plotly_white",
    )

    return fig


# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container(
    [
        html.H1("PCR SIM DASHBOARD"),
        html.H5("by K. Soikum, S.Soikum and L.Thomsen, 2024"),
        dbc.Row(
            [
                dcc.Dropdown(
                    id="fasta-dropdown",
                    options=[
                        {"label": "Ebola Virus (Zaire ebolavirus): NC_002549", "value": "NC_002549"},
                        {"label": "Marburg Virus (Marburg marburgvirus): NC_001608", "value": "NC_001608"},
                        {"label": "H5N1 Avian Influenza (Influenza A virus): CY115771", "value": "CY115771"},
                        {"label": "HIV-1 (Human Immunodeficiency Virus): NC_001802", "value": "NC_001802"},
                        {"label": "Human telomerase reverse transcriptase: NM_198253", "value": "NM_198253"},
                        {"label": "Bioluminescent enzyme (Luciferase) from firefly: NM_001293166", "value": "NM_001293166"},
                        {"label": "Human programmed cell death protein 1 (PD-1): NM_005018", "value": "NM_005018"},
                    ],
                    placeholder="Enter or select a FASTA ID...",
                    className="mb-3",
                    searchable=True,
                    clearable=True,
                    style={"width": "100%"},
                ),
            ],
            className="mb-3",
        ),
        dbc.Button("Fetch and Design Primers", id="fetch-button", color="primary", className="mb-3"),
        dbc.Tabs(
            [
                dbc.Tab(label="Visualization", tab_id="visualization"),
                dbc.Tab(label="Primer Table", tab_id="primer-table"),
                dbc.Tab(label="PCR Simulations", tab_id="pcr-simulations"),
                dbc.Tab(label="Capillary Electrophoresis", tab_id="capillary-electrophoresis"),
            ],
            id="tabs",
            active_tab="visualization",
        ),
        dbc.Spinner(
            [
                dcc.Store(id="store"),
                html.Div(id="tab-content", className="p-4"),
            ],
            delay_show=100,
        ),
    ],
    fluid=True,
)



#@app.callback(
#    Output("fasta-id", "value"),
#    [Input("dropdown", "value")]
#)
#def populate_input_from_dropdown(selected_value):
#    if selected_value:
#        return selected_value
#    return ""


@app.callback(
    Output("tab-content", "children"),
    [Input("tabs", "active_tab"), Input("store", "data")]
)
def render_tab_content(active_tab, data):
    """
    Render the content dynamically for the active tab.
    """
    if data is None:
        return "No content available. Please fetch data by entering a valid FASTA ID."

    if active_tab:
        if active_tab == "visualization":
            figure = data.get("visualization", None)
            if not figure:
                return "No visualization available."
            return dcc.Graph(figure=figure)
        elif active_tab == "primer-table":
            table_data = data.get("table", [])
            if not table_data:
                return "No primer data available."
            return dash_table.DataTable(
                id="primer-table",
                data=table_data,
                columns=[
                    {"name": i, "id": i}
                    for i in (table_data[0].keys() if table_data else [])
                ],
                style_table={"overflowX": "auto"},
            )
        elif active_tab == "pcr-simulations":
            amplification = data["pcr"]["amplification"]
            melting = data["pcr"]["melting"]

            amplification_fig = go.Figure()
            for curve in amplification:
                amplification_fig.add_trace(go.Scatter(x=curve["x"], y=curve["y"], name=curve["name"]))

            melting_fig = go.Figure()
            for curve in melting:
                melting_fig.add_trace(go.Scatter(x=curve["x"], y=curve["y"], name=curve["name"]))

            return html.Div([
                dcc.Graph(figure=amplification_fig),
                dcc.Graph(figure=melting_fig)
            ])
        elif active_tab == "capillary-electrophoresis":
            # Fetch data for Capillary Electrophoresis
            ce_figure = data["capillary-electrophoresis"] if data and "capillary-electrophoresis" in data else go.Figure()
            return dcc.Graph(figure=ce_figure)

    return "No tab selected."



@app.callback(
    Output("store", "data"),
    [Input("fetch-button", "n_clicks")],
    [State("fasta-dropdown", "value")]
)
def generate_data(n_clicks, fasta_id):
    if not n_clicks or not fasta_id:
        return {"visualization": go.Figure(), "table": [], "pcr": {"amplification": [], "melting": []}}

    try:
        sequence = fetch_sequence_by_fasta(fasta_id)
        
        # Fetch primers and full results from Primer3
        full_results = primer3.bindings.design_primers(
            {
                'SEQUENCE_ID': 'example',
                'SEQUENCE_TEMPLATE': sequence,
                'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],
            },
            {
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
        )

        # Extract primers in a user-friendly format
        primers = [
            {
                "Index": i + 1,
                "Left Primer": full_results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                "Right Primer": full_results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                "Penalty": full_results[f'PRIMER_PAIR_{i}_PENALTY'],
                "Left Position": full_results[f'PRIMER_LEFT_{i}'][0],
                "Right Position": full_results[f'PRIMER_RIGHT_{i}'][0],
                "Product Size": full_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            }
            for i in range(full_results['PRIMER_LEFT_NUM_RETURNED'])
        ]

        # Generate visualization
        fig = plot_sequence_with_primers(sequence, primers, max_sequence_length=50000)


        # Pass full_results to simulate_pcr_amplification
        amplification_curves, melting_curves = simulate_pcr_amplification(primers, sequence, full_results)
        # Simulate Capillary Electrophoresis
        amplicon_lengths = [primer["Product Size"] for primer in primers]
        ce_figure = simulate_capillary_electrophoresis(amplicon_lengths)
        return {
            "visualization": fig.to_dict(),
            "table": primers,
            "pcr": {"amplification": amplification_curves, "melting": melting_curves},
            "capillary-electrophoresis": ce_figure.to_dict(),  # Add CE data
        }
    except Exception as e:
        return {"visualization": go.Figure().to_dict(), "table": [{"Error": str(e)}], "pcr": {"amplification": [], "melting": []}}



# Run the server
if __name__ == "__main__":
    app.run_server(debug=True)
