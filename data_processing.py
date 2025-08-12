"""
Data processing utilities for gene expression analysis.
Contains functions for filtering, aggregating, and processing gene expression data from both WormSeq and CenGen datasets.
"""

import pandas as pd
import os
from pathlib import Path

#Configuration constants
SCALED_TPM_THRESHOLD = 400
SPERMATHECA_REGIONS = {
    "Spermatheca-Uterine junction" : "valve",
    "Spermatheca neck distal": "neck", 
    "Spermatheca neck proximal": "neck", 
    "Spermatheca bag distal": "bag", 
    "Spermatheca bag proximal": "bag"
}

def filter_high_expression_genes(input_folder, expression_threshold=SCALED_TPM_THRESHOLD, expression_column="scaled_TPM"):
    """
    Filter genes with expression above threshold and sort by expression level
    """
    for filename in os.listdir(input_folder):
        if filename.endswith("csv") and not filename.endswith("_filtered.csv"):
            input_path = Path(input_folder) / filename
            df = pd.read_csv(input_path)

            filtered_df = df[df[expression_column] > expression_threshold].copy()
            filtered_df.sort_values(by=expression_column, ascending=False, inplace=True)

            output_filename = filename.replace(".csv", "_filtered.csv")
            output_path = input_path.parent / output_filename
            filtered_df.to_csv(output_path, index=False)

            print(f"{output_filename}: Filtered {len(filtered_df)} genes with {expression_column} > {expression_threshold}")

def aggregate_expression_matrix(input_folder, output_folder, gene_column, expression_column, output_filename):
    """
    Aggregate expression data into a gene x region matrix
    """

    os.makedirs(output_folder, exist_ok=True)
    
    region_expression_data = {}

    #process each filtered csv file
    for filename in os.listdir(input_folder):
        if filename.endswith("_filtered.csv"):
            region_name = filename.replace("_filtered.csv", "")
            filepath = Path(input_folder) / filename
            df = pd.read_csv(filepath)

            region_expression_data[region_name] = df.set_index(gene_column)[expression_column]

    expression_matrix = pd.DataFrame(region_expression_data)
    expression_matrix.dropna(how='all', inplace=True)

    output_path = Path(output_folder) / output_filename
    expression_matrix.to_csv(output_path)

    print(f"Saved expression matrix to: {output_path}")
    
def process_wormseq_data():
    """
    Process WormSeq spermatheca data through filtering and aggregation
    """
    filter_high_expression_genes("spermatheca_cell_types")
    aggregate_expression_matrix(
        input_folder="spermatheca_cell_types",
        output_folder="scaled_TPM_heatmap",
        gene_column="gene_short_name",
        expression_column="scaled_TPM", 
        output_filename="scaled_TPM_matrix.csv"
    )

def process_cengen_data():
    """
    Process CenGen data through filtering and aggregation
    """
    filter_high_expression_genes("cengen", expression_column="Expression level")
    aggregate_expression_matrix(
        input_folder="cengen",
        output_folder="scaled_TPM_heatmap", 
        gene_column="Gene name", 
        expression_column="Expression level", 
        output_filename="cengen_scaled_TPM_matrix.csv"
    )

if __name__ == "__main__":
    process_wormseq_data()
    process_cengen_data()