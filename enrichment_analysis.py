"""
Filter and prepare gene lists for GO enrichment analysis.
Filteres gene based on spermtheca-specific expression patterns and prepares them for downstream ontological enrichment analysis. 
"""

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


SPERMATHECA_COMPONENTS = {
    "bag": ["Spermatheca bag distal", "Spermatheca bag proximal"], 
    "neck": ["Spermatheca neck distal", "Spermatheca neck proximal"], 
    "valve": ["Spermatheca-Uterine junction"]
}

def determine_spermatheca_component(filename):
    """
    Determine which spermatheca component a file represents based on filename
    """
    filename_lower = filename.lower()
    for component, region_keywords in SPERMATHECA_COMPONENTS.items():
        if any(keyword.lower() in filename_lower for keyword in region_keywords):
            return component
    return None

def filter_spermatheca_specific_genes(input_folder="w_relative_TPM", output_folder="filtered_for_enrichment"):
    """
    Filter gene for enrichment analysis based on spermatheca specific expression. 
    """
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if not filename.endswith(".csv"):
            continue
        component = determine_spermatheca_component(filename)
        if component is None:
            print(f"Skipping {filename} - no component mapping found")
            continue
        region_name = filename.replace(".csv", "")
        input_path = Path(input_folder) / filename
        df = pd.read_csv(input_path)

        spermatheca_specific = df[df['max_in_spermatheca'] == True].copy()
        spermatheca_output = Path(output_folder) / f"{region_name}_filtered_spermatheca.csv"
        spermatheca_specific.to_csv(spermatheca_output, index=False)

        component_specific = spermatheca_specific[
            spermatheca_specific['max_in_same_component'] == True
        ].copy()
        component_output = Path(output_folder) / f"{region_name}_filtered_component.csv"
        component_specific.to_csv(component_output, index=False)
        
        print(f"Processed {region_name}: ")
        print(f"Spermatheca-specific genes: {len(spermatheca_specific)}")
        print(f" Component-specific genes: {len(component_specific)}")

def combine_enrichment_results(enrichment_folder="enrichment_results", output_file="combined_spermatheca_enrichment.csv"):
    """
    Combine GO enrichment results from multiple spermatheca regions.
    """
    spermatheca_files = [
         f for f in os.listdir(enrichment_folder)
         if f.endswith("_enrichment_spermatheca.csv")
    ]
    if not spermatheca_files:
        print("No spermatheca enrichment files found")
        return
    def extract_region_name(filename):
        """
        Extract region name from enrichment filename.
        """
        return filename.replace("_enrichment_spermatheca.csv", "")
    combined_results = []

    for filename in spermatheca_files:
        filepath = Path(enrichment_folder) / filename
        region_name = extract_region_name(filename)

        enrichment_df = pd.read_csv(filepath)
        enrichment_df["region"] = region_name
        combined_results.append(enrichment_df)

    final_combined = pd.concat(combined_results, ignore_index=True)
    final_combined.to_csv(output_file, index=False)

    print(f"Combined enrichment results saved to: {output_file}")
    print(f"Regoins included: {final_combined['region'].unique()}")
    print(f"Total GO terms: {final_combined['name'].nunique()}")
    return final_combined
    
def create_relative_tpm_heatmap(input_folder="w_relative_TPM", output_folder="TPM_heatmap"):
    """
    Create heatmap visualization of relative TPM values across spermatheca regions.
    """
    os.makedirs(output_folder, exist_ok=True)
    region_label_mapping = {
        "Spermatheca-Uterine junction": "valve", 
        "Spermatheca neck distal": "neck_distal", 
        "Spermatheca neck proximal": "neck_proximal", 
        "Spermatheca bag proximal": "bag_proximal", 
        "Spermatheca bag distal": "bag_distal"
    }
    region_dataframes = []

    for filename in os.listdir(input_folder):
        if not filename.endswith(".csv"):
            continue
        region_full_name = filename.replace(".csv", "")
        region_short_name = region_label_mapping.get(region_full_name)
        if not region_short_name:
            print(f"Skipping unknown region: {filename}")
            continue
        filepath = Path(input_folder) / filename
        df = pd.read_csv(filepath)

        spermatheca_max_genes = df[df['max_in_spermatheca'] == True].copy()

        processed_df = spermatheca_max_genes[["gene_ID", "gene_short_name", "relative_TPM"]].copy()
        processed_df.rename(columns={"relative_TPM": region_short_name}, inplace=True)
        region_dataframes.append(processed_df)

    if not region_dataframes:
        print("No valid region data found")
        return
    merged_matrix = region_dataframes[0]
    for df in region_dataframes[1:]:
        merged_matrix = pd.merge(
            merged_matrix, df, 
            on=["gene_ID", "gene_short_name"], 
            how="outer"
        )
    matrix_output = Path(output_folder) / "relative_TPM_matrix.csv"
    merged_matrix.to_csv(matrix_output, index=False)
    heatmap_data = merged_matrix.set_index("gene_short_name")
    expression_columns = [col for col in heatmap_data.columns if col != "gene_ID"]
    plot_data = heatmap_data[expression_columns]

    plt.figure(figsize=(10, max(8, 0.25 * len(plot_data))))
    sns.heatmap(
        plot_data,
        cmap="YlGnBu", 
        linewidths=0.5, 
        linecolor='gray', 
        cbar_kws={'label': 'Relative TPM'}, 
        xticklabels=True, 
        yticklabels=True
    )
    plt.title("Relative TPM Across Spermatheca Regions")
    plt.xlabel("Spermatheca Region")
    plt.ylabel("Gene")
    plt.tight_layout()

    heatmap_output = Path(output_folder) / "relative_TPM_heatmap.png"
    plt.savefig(heatmap_output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Relative TPM analysis complete. Results saved to: {output_folder}")

def main():
    """
    Main pipeline for enrichment filtering and analysis.
    """
    print("Filtering genes for enrichment analysis...")
    filter_spermatheca_specific_genes()
    print("Creating relative TPM heatmap...")
    create_relative_tpm_heatmap()
    print("Enrichment analysis pipeline complete!")

if __name__ == "__main__":
    main()