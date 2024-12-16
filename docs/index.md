---
​---
title: Home
permalink: /PCRSIM/
​---
---

# Title: PCR SIM DOCUMENTATION

#### Authors: Kawin Soikum, Soiwisa Soikum, and Lars Thomsen 

Affiliation: [GNACODE.INC](www.gnacode.com), Medicine Hat, Alberta, Canada

Email for correspondence: [**Lars Thomsen**](mailto:lt@gnacode.com)

This document provides a detailed explanation of the code's functionality, which implements a dashboard for PCR simulation and DNA primer design. It is structured to help developers and users understand the purpose and operation of each function.

![](.\images\sequence_view.png)

## Overview

The code:
1. Fetches DNA sequences from NCBI using a FASTA ID.
2. Designs primers for PCR using Primer3.
3. Simulates PCR amplification and melting curves.
4. Visualizes DNA sequences and primer locations.
5. Simulates capillary electrophoresis for amplicons.

### Libraries Used
- `os`, `sys`: Handle file system and script paths.
- `primer3`: Designs primers for DNA sequences.
- `dash`, `dash_table`, `dcc`, `html`: Build an interactive web-based dashboard.
- `dash_bootstrap_components`: Adds Bootstrap styling to the dashboard.
- `Bio` (`Entrez`, `SeqIO`, `Seq`): Fetches and processes DNA sequences.
- `plotly.graph_objects`: Creates interactive plots.
- `matplotlib.colors`, `itertools`, `numpy`: Perform data visualization and mathematical operations.

---

## Functions

### 1. **`clean_dna_sequence(dna_sequence)`**
Cleans and validates a DNA sequence by removing whitespace and ensuring it contains only valid bases (A, T, C, G).

#### Args:
- `dna_sequence` (str): Raw DNA sequence.

#### Returns:
- Cleaned DNA sequence (str).

#### Raises:
- `ValueError`: If the sequence contains invalid characters.

---

### 2. **`fetch_sequence_by_fasta(fasta_id)`**
Fetches a DNA sequence from the NCBI nucleotide database using a FASTA ID.

#### Args:
- `fasta_id` (str): NCBI FASTA ID.

#### Returns:
- DNA sequence (str).

#### Notes:
- Requires the user to provide a valid email address (`Entrez.email`).

---

### 3. **`design_primers(sequence)`**
Designs PCR primers for a given DNA sequence using Primer3.

#### Args:
- `sequence` (str): DNA sequence.

#### Returns:
- List of dictionaries containing primer details (left primer, right primer, penalties, positions, product size).

---

### 4. **`add_product_region_shapes(shapes, product_start, product_end, bases_per_row, light_color)`**
Adds rectangle shapes representing product regions to a DNA visualization.

#### Args:
- `shapes` (list): List of shape dictionaries for Plotly.
- `product_start`, `product_end` (int): Start and end positions of the product region.
- `bases_per_row` (int): Number of bases displayed per row.
- `light_color` (str): Color for the product region.

---

### 5. **`plot_sequence_with_primers(sequence, primer_pairs, target_aspect_ratio=1.5, max_sequence_length=50000)`**
Visualizes a DNA sequence and highlights primer positions using Plotly.

#### Args:
- `sequence` (str): DNA sequence.
- `primer_pairs` (list): List of primers.
- `target_aspect_ratio` (float): Aspect ratio for visualization.
- `max_sequence_length` (int): Maximum sequence length for full visualization.

#### Returns:
- Plotly figure object.

---

### 6. **`reverse_complement(sequence)`**
Generates the reverse complement of a DNA sequence.

#### Args:
- `sequence` (str): DNA sequence.

#### Returns:
- Reverse complement of the DNA sequence (str).

---

### 7. **`calculate_tm(sequence)`**
Calculates the melting temperature (Tm) of a DNA sequence using the nearest-neighbor model.

#### Args:
- `sequence` (str): DNA sequence.

#### Returns:
- Melting temperature (float).

---

### 8. **`adjust_efficiency(original_efficiency, product_length)`**
Adjusts the amplification efficiency based on the product length.

#### Args:
- `original_efficiency` (float): Initial efficiency.
- `product_length` (int): Length of the PCR product.

#### Returns:
- Adjusted efficiency (float).

---

### 9. **`find_product_length(primer_pair, sequence)`**
Calculates the length of the PCR product based on primer positions.

#### Args:
- `primer_pair` (tuple): Pair of primers (forward, reverse).
- `sequence` (str): DNA sequence.

#### Returns:
- Product length (int).

#### Raises:
- `ValueError`: If primers are not found in the sequence.

---

### 10. **`simulate_pcr_amplification(primer_pairs, sequence, results)`**
Simulates PCR amplification and melting curves for given primers and sequences.

#### Args:
- `primer_pairs` (list): List of primer pairs.
- `sequence` (str): DNA sequence.
- `results` (dict): Primer3 results.

#### Returns:
- Amplification curves (list).
- Melting curves (list).

---

### 11. **`simulate_capillary_electrophoresis(amplicon_lengths)`**
Simulates capillary electrophoresis for given amplicon lengths.

#### Args:
- `amplicon_lengths` (list): List of amplicon lengths.

#### Returns:
- Plotly figure object for electrophoresis.

---

## Dash Application

### Layout
The dashboard includes:
1. **Dropdown Menu**: Select or enter a FASTA ID.
2. **Tabs**: 
   - Visualization: DNA and primers.
   - Primer Table: Primer data.
   - PCR Simulations: Amplification and melting curves.
   - Capillary Electrophoresis: Electrophoresis plot.

### Callbacks
1. **Fetch Data**: Fetches sequence and designs primers on button click.
2. **Render Tab Content**: Dynamically renders content for the selected tab.
3. **Store Data**: Saves visualization, primer data, and simulations for rendering.

---

## Execution

 The programs are written in Python 3.11. and the requirements.txt file contains a list over the Python libraries and their versions required to run the program. If the programs are downloaded and installed using Git (see Installation below)

[**INSTALLATION**](Installation.md)

[**LIVESHARE WITH VSCODE AND GITHUB CODESPACES**](LiveShare.md)

[**MIT LICENSE**](MIT-license.md) for free use

[**ACKNOWLEDGEMENT**](Acknowledgement.md) on use





