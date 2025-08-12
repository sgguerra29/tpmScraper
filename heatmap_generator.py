"""
Generate heatmaps for gene expression data visualization.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

def create_expression_heatmap(input_file, output_folder, title, output_filename, expression_column="scaled_TPM", log_transform=True):
    """
    Create a heatmap from gene expression matrix data.
    """
    expression_df = pd.read_csv(input_file, index_col=0)

    if log_transform:
        plot_data = np.log10(expression_df + 1)
        colorbar_label = 'log10(scaled_TPM +1)'
    else:
        plot_data = expression_df
        colorbar_label = 'scaled_TPM'

    height = min(0.2 * len(plot_data), 60)

    plt.figure(figsize=(10, height))
    sns.heatmap(
        plot_data, 
        cmap="viridis", 
        linewidths=0.5, 
        linecolor='gray', 
        cbar_kws={'label': colorbar_label}, 
        xticklabels=True, 
        yticklabels=True
    )

    plt.title(title, fontsize=14)
    plt.xlabel("Region")
    plt.ylabel("Gene")
    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    output_path = Path(output_folder) / output_filename
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved heatmap tp: {output_path}")

def generate_wormseq_heatmap():
    """
    Generate heatmap for WormSeq spermatheca expression data
    """
    create_expression_heatmap(
        input_file="scaled_TPM_heatmap/scaled_TPM_matrix.csv", 
        output_folder="scaled_TPM_heatmap", 
        title="Gene Expression Across Spermatheca Regions", 
        output_filename="wormseq_scaled_TPM_heatmap.png"
    )

def generate_cengen_heatmap():
    """
    Generate heatmap for CenGen expression data
    """
    create_expression_heatmap(
        input_file="scaled_TPM_heatmap/cengen_scaled_TPM_matrix.csv", 
        output_folder="scaled_TPM_heatmap", 
        title="Gene Expression in CenGen Regions", 
        output_filename="cengen_scaled_TPM_heatmap.png"
    )

if __name__ == "__main__":
    generate_cengen_heatmap()
    generate_wormseq_heatmap()
