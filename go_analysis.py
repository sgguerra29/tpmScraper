"""
Gene Ontology enrichment analysis and visualization. 
Performs enrichment analysis using g:Profiler and creates plots to visualize GO term distributions across spermatheca regions. 
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from gprofiler import GProfiler
from pathlib import Path

SPERMATHECA_REGIONS = [
    "neck distal", "neck proximal", 
    "bag distal", "bag proximal", 
    "Spermatheca-Uterine junction"
]

REGION_GROUPS = {
    "neck": ["neck distal", "neck proximal"], 
    "bag": ["bag distal", "bag proximal"], 
    "valve": ["Spermatheca-Uterine junction"]
}

def perform_go_enrichment_analysis(input_folder = "filtered_for_enrichment", output_folder="enrichment_results"):
    """
    Perform GO enrichment analysis on filtered gene sets
    """
    os.makedirs(output_folder, exist_ok=True)

    profiler = GProfiler(return_dataframe=True)

    for filename in os.listdir(input_folder):
        if filename.endswith("_filtered_spermatheca.csv"):
            region_name = filename.replace("_filtered_spermatheca.csv", "")
            process_enrichment_file(filename, region_name, "spermatheca", input_folder, output_folder, profiler)
        elif filename.endswith("_filtered_component.csv"):
            region_name = filename.replace("_filtered_component.csv", "")
            process_enrichment_file(filename, region_name, "component", input_folder, output_folder, profiler)


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
    print(f"Regions included: {final_combined['region'].unique()}")
    print(f"Total GO terms: {final_combined['name'].nunique()}")
    return final_combined

def process_enrichment_file(filename, region_name, analysis_type, input_folder, output_folder, profiler):
    """
    Process a single gene list file for GO enrichment
    """
    input_path = Path(input_folder) / filename
    print(f"Processing {region_name} - {analysis_type} analysis")
   
    df = pd.read_csv(input_path)
    gene_list = df['gene_short_name'].dropna().unique().tolist()
    if not gene_list:
        print(f"No genes found in {filename}")
        return
    try:
        enrichment_results = profiler.profile(
            organism='celegans',
            query=gene_list, 
            user_threshold=0.05, 
            sources=["GO:BP", "GO:MF", "GO:CC"], 
            no_evidences=False
        )
        enrichment_results['region'] = region_name
        output_filename = f"{region_name}_enrichment_{analysis_type}.csv"
        output_path = Path(output_folder) / output_filename
        enrichment_results.to_csv(output_path, index=False)

        print(f"Saved {len(enrichment_results)} GO terms for {region_name}")
    except Exception as e:
        print(f"Error processing {region_name}: {e}")

def main():
    """
    Main analysis pipeline
    """
    perform_go_enrichment_analysis()
    print("Combining enrichment results...")
    combine_enrichment_results()
    print("GO enrichment analysis complete - gene list ready for extraction")

if __name__ == "__main__":
    main()





