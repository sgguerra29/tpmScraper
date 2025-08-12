"""
Extract gene lists from GO term enrichment resuls by spermatheca region.
This script processes combined GO enrichment data and creates separate CSV files for each spermatheca region containg GO terms and their associated genes. 
"""

import pandas as pd
from pathlib import Path

def clean_go_id(go_id):
    """
    Clean GO term ids for use in filename by replacing colons with underscores:
    """
    return go_id.replace(":", "_")

def clean_region_name(region_name):
    """
    Clean region names to match expected format for plotting. 
    """
    if region_name.startswith("Spermatheca "):
        return region_name.replace("Spermatheca ", "")
    return region_name

def extract_gene_list_by_region(input_file="combined_spermatheca_enrichment.csv", output_folder="go_genes_list"):
    """
    Extract and save gene lists grouped by spermatheca region.
    """
    output_path = Path(output_folder)
    output_path.mkdir(exist_ok=True)

    enrichment_df = pd.read_csv(input_file)

    #Process each region separately
    for region in enrichment_df["region"].unique():
        region_data = enrichment_df[enrichment_df["region"] == region]
        gene_data_list = []
        for _, row in region_data.iterrows():
            go_id = row["native"]
            description = row["name"]
            genes = row["intersections"].split(",") if pd.notna(row["intersections"]) else []
        
            gene_data_list.append({
                "GO_ID" : go_id, 
                "Description" : description, 
                "Gene_Count" : len(genes), 
                "Genes": ";".join(genes)
        })

        cleaned_region = clean_region_name(region)
        region_output_path = output_path / f"{cleaned_region}_go_genes_lists.csv"
        region_df = pd.DataFrame(gene_data_list)
        region_df.to_csv(region_output_path, index=False)

    print(f"Go term gene lists extracted for {len(enrichment_df['region'].unique())} regions")


if __name__ == "__main__":
    extract_gene_list_by_region()


