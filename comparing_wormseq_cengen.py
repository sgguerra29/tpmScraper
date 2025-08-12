"""
Compare and visualize gene expression data between WormSeq and CenGen datasets. 
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from matplotlib_venn import venn2

def normalize_dataset_columns(df, source_name):
    """
    Normalize column names across different datasets.
    """
    df = df.copy()
    
    if "cengen" in source_name.lower():
        df = df.rename(columns={
            "Gene name": "gene",
            "Expression level": "scaled_TPM"
        })
    elif source_name == "wormseq spermatheca":
        df = df.rename(columns={
            "gene_short_name": "gene", 
            "max_scaled_TPM": "scaled_TPM"
        })
    elif source_name == "wormseq sp_ut":
        df = df.rename(columns={"gene_short_name": "gene"})
        df = df.groupby("gene", as_index=False)["scaled_TPM"].max()

    df = df[["gene", "scaled_TPM"]].copy()
    df["source"] = source_name

    return df

def combine_expression_datasets():
    """
    Combine WormSeq and CenGEn expression data into a single comparison table.
    """
    dataset_files = {
        "cengen sp_ut": "cengen/CenGen sp_ut_filtered.csv", 
        "cengen spermatheca": "cengen/CenGen spermatheca_filtered.csv", 
        "wormseq sp_ut": "spermatheca_cell_types/Spermatheca-Uterine junction_filtered.csv", 
        "wormseq spermatheca": "spermatheca_cell_types/merged_wormseq_spermatheca.csv"
    }
    combined_data = []
    for source_name, file_path in dataset_files.items():
        if not Path(file_path).exists():
            print(f"Warning: File not found - {file_path}")
            continue
        df = pd.read_csv(file_path)
        normalized_df = normalize_dataset_columns(df, source_name)
        combined_data.append(normalized_df)

    all_expression_data = pd.concat(combined_data, ignore_index=True)
    comparison_matrix = all_expression_data.pivot_table(
        index="gene", 
        columns="source", 
        values="scaled_TPM"
    )

    output_path = "gene_expression_comparison.csv"
    comparison_matrix.to_csv(output_path)
    print(f"Dataset comparison saved to: {output_path}")
    print(f"Total genes in comparison: {len(comparison_matrix)}")
    return comparison_matrix

def create_dataset_comparison_plots(comparison_df, output_folder="combined_datasets"):
    """
    Generate visualization plots comparing datasets.
    """
    os.makedirs(output_folder, exist_ok=True)
    spermatheca_common = comparison_df.dropna(subset=["cengen spermatheca", "wormseq spermatheca"])
    sput_common = comparison_df.dropna(subset=["cengen sp_ut", "wormseq sp_ut"])
    common_genes_df = pd.concat([spermatheca_common, sput_common]).drop_duplicates()
    common_genes_df.to_csv("common_genes_by_regions.csv", index=False)
    create_common_genes_heatmap(common_genes_df, output_folder)
    create_dataset_venn_diagram(comparison_df, output_folder)
    create_expression_scatter_plot(common_genes_df, output_folder)

def create_common_genes_heatmap(common_genes_df, output_folder):
    """
    Create heatmap of expression levels for genes common to both datasets.
    """
    heatmap_data = common_genes_df.fillna(0)
    log_transformed = np.log10(heatmap_data + 1)
    
    height = min(0.2 * len(log_transformed), 60)
    plt.figure(figsize=(10, height))
    sns.heatmap(
        log_transformed, 
        cmap="viridis", 
        cbar_kws={'label': 'log10(scaled_TPM + 1)'}, 
        linewidths=0.1, 
        linecolor='white'
    )
    plt.title("Expression of Common Genes Across Datasets", fontsize=14)
    plt.xlabel("Dataset & Region")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "heatmap_common_genes.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_dataset_venn_diagram(comparison_df, output_folder):
    """
    Create Venn diagram showing gene overlap between datasets.
    """
    cengen_expressed = set(comparison_df[
        (~comparison_df["cengen sp_ut"].isna()) |
        (~comparison_df["cengen spermatheca"].isna())].index)
    wormseq_expressed = set(comparison_df[
        (~comparison_df["wormseq sp_ut"].isna()) |
        (~comparison_df["wormseq spermatheca"].isna())].index)

    plt.figure(figsize=(8, 6))
    venn2([cengen_expressed, wormseq_expressed], set_labels=("CenGen", "WormSeq"))
    plt.title("Gene Coverage Overlap Between Datasets")
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "venn_diagram_overlap.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_expression_scatter_plot(common_genes_df, output_folder):
    """
    Create scatter plot showing expression patterns across datasets.
    """
    plot_data = common_genes_df.reset_index().melt(
        id_vars="gene", 
        var_name="dataset_region", 
        value_name="scaled_TPM"
    )
    plot_data = plot_data.dropna(subset=["scaled_TPM"])
    plt.figure(figsize=(12, max(8, 0.3 * len(common_genes_df))))
    sns.scatterplot(
        data=plot_data, 
        x="dataset_region", 
        y="gene",
        size="scaled_TPM", 
        hue="scaled_TPM", 
        palette="viridis", 
        sizes=(20, 200), 
        alpha=0.7 
    )

    plt.title("Gene Expression Levels by Dataset & Region")
    plt.xlabel("Dataset & Region")
    plt.ylabel("Gene")
    plt.xticks(rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "expression_scatter_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """
    Main comparison analysis pipeline.
    """
    print("Combining expression datasets...")
    comparison_matrix = combine_expression_datasets()
    print("Creating comparison visualizations....")
    create_dataset_comparison_plots(comparison_matrix)
    print("Dataset comparison analysis complete!")

if __name__ == "__main__":
    main()