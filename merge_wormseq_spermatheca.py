"""
Merge WormSeq data across spermatheca regions and find max scaled_TPM per gene.
Combines expression data from multiple spermatheca subregions and identifies max scaled_TPM for each gene across all regions for downstream comparisons with CenGen data
"""


import os
import pandas as pd
from pathlib import Path

def merge_spermatheca_regions(input_folder="spermatheca_cell_types", output_filename="merged_wormseq_spermatheca.csv"):
    """
    Merge scaled_TPM values across spermatheca regions, keeping max per gene. 
    """

    region_files = [
    "Spermatheca neck distal_filtered.csv", 
    "Spermatheca neck proximal_filtered.csv",
    "Spermatheca bag distal_filtered.csv",
    "Spermatheca bag proximal_filtered.csv"
    ]

    gene_max_expression = {}

    #Process each region file
    for filename in region_files:
        filepath = Path(input_folder) / filename
        if not filepath.exists():
            print (f"Warning: File not found - {filepath}")
            continue

        region_df = pd.read_csv(filepath)

        #update max scaled TPM for each gene
        for _, row in region_df.iterrows():
            gene_name = row["gene_short_name"]
            scaled_TPM = row["scaled_TPM"]

            if (gene_name not in gene_max_expression or scaled_TPM > gene_max_expression[gene_name]):
                gene_max_expression[gene_name] = scaled_TPM

    merged_data = pd.DataFrame({
        "gene_short_name" : list(gene_max_expression.keys()),
        "max_scaled_TPM" : list(gene_max_expression.values())
    })

    merged_data = merged_data.sort_values("max_scaled_TPM", ascending=False)

    output_path = Path(input_folder) / output_filename
    merged_data.to_csv(output_path, index=False)

    print(f"Merged data saved to: {output_path}")
    print(f"Total genes processed: {len(merged_data)}")
    print(f"\nTop 5 highly expressed genes:")
    print(merged_data.head())

    return merged_data

if __name__ == "__main__":
    merge_spermatheca_regions()

