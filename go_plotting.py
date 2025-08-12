"""
Create GO term visualization plots from extracted gene lists. 
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
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

def load_go_gene_data(input_folder):
    """
    Load GO term data from region-specific CSV files.
    """
    region_go_mapping = {}
    for region in SPERMATHECA_REGIONS:
        file_path = Path(input_folder) / f"{region}_go_genes_lists.csv"
        if not file_path.exists():
            print(f"File not found - {file_path}")
            continue
        df = pd.read_csv(file_path)
        go_term_genes = {}

        for _, row in df.iterrows():
            go_term = row['Description']
            genes = row['Genes'].split(';') if pd.notna(row['Genes']) else []
            go_term_genes[go_term] = genes
        region_go_mapping[region] = go_term_genes
    return region_go_mapping

def compute_go_term_stats(region_go_mapping):
    """
    Compute gene counts and unique GO terms per region. 
    """
    region_gene_counts = {}
    for region, go_mapping in region_go_mapping.items():
        go_gene_counts = {go_term: len(genes) for go_term, genes in go_mapping.items()}
        region_gene_counts[region] = go_gene_counts
    all_go_terms = defaultdict(set)
    for region, go_mapping in region_go_mapping.items():
        for go_term in go_mapping:
            all_go_terms[go_term].add(region)
    region_unique_terms = defaultdict(list)
    for go_term, regions in all_go_terms.items():
        if len(regions) == 1:
            region = next(iter(regions))
            region_unique_terms[region].append(go_term)
    return region_gene_counts, region_unique_terms

def create_go_term_plots(region_go_mapping, output_folder="go_term_plots"):
    """
    Generate comprehensive GO term visualization plots. 
    """
    os.makedirs(output_folder, exist_ok=True)
    region_gene_counts, region_unique_terms = compute_go_term_stats(region_go_mapping)
    plot_top_go_terms_per_region(region_gene_counts, output_folder)
    plot_unique_go_terms_comparison(region_unique_terms, output_folder)
    create_stacked_comparison_plots(region_go_mapping, output_folder)
    plot_grouped_tissue_analysis(region_go_mapping, output_folder)

def plot_top_go_terms_per_region(region_gene_counts, output_folder, top_n=100):
    """
    Create bar plots for top GO terms by gene count for each region.
    """
    for region, go_counts in region_gene_counts.items():
        if not go_counts:
            continue

        sorted_terms = sorted(go_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
        go_terms, gene_counts = zip(*sorted_terms)

        plt.figure(figsize=(12, 6))
        plt.bar(range(len(go_terms)), gene_counts, color='steelblue')
        plt.title(f"Top {top_n} GO Terms by Gene Count - {region}")
        plt.ylabel("Number of Genes")
        plt.xlabel("GO Terms")
        plt.xticks(range(len(go_terms)), go_terms, rotation=45, ha='right')
        plt.tight_layout()

        safe_region_name = region.replace(' ', '_').replace('-', '_')
        plt.savefig(Path(output_folder) / f"{safe_region_name}_top_go_terms.png", dpi=300, bbox_inches='tight')
        plt.close()

def plot_unique_go_terms_comparison(region_unique_terms, output_folder):
    """
    Create bar chart comparing number of unique GO terms per region.
    """
    unique_counts = [len(region_unique_terms.get(region, [])) for region in SPERMATHECA_REGIONS]
    plt.figure(figsize=(10, 6))
    bars = plt.bar(SPERMATHECA_REGIONS, unique_counts, color='coral')
    plt.title("Unique GO Terms per Spermatheca Region")
    plt.ylabel("Number of Unique GO Terms")
    plt.xlabel("Region")
    plt.xticks(rotation=45, ha='right')

    for bar, count in zip(bars, unique_counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, str(count), ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "unique_go_terms_per_region.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_stacked_comparison_plots(region_go_mapping, output_folder):
    """
    Create stacked bar charts for specific region comparisons.
    """
    comparisons = [
        (["neck distal", "neck proximal"], "neck_distal_vs_proximal.png"), 
        (["bag distal", "bag proximal"], "bag_distal_vs_proximal.png"), 
        (SPERMATHECA_REGIONS, "all_regions_stacked.png")
    ]
    for regions, filename in comparisons:
        create_stacked_bar_plot(region_go_mapping, output_folder, regions, filename)

def create_stacked_bar_plot(region_go_mapping, output_folder, region_list, filename):
    """
    Create a stacked bar chart for specified regions. 
    """
    go_term_counts = defaultdict(dict)

    for region in region_list:
        if region not in region_go_mapping:
            continue
        for go_term, genes in region_go_mapping[region].items():
            go_term_counts[go_term][region] = len(genes)

    df = pd.DataFrame(go_term_counts).fillna(0).T
    df = df.reindex(columns=region_list)

    plt.figure(figsize=(14, 8))
    df.plot(kind='bar', stacked=True, ax=plt.gca())
    plt.title(f"GO Term Gene Counts: {' vs '.join(region_list)}")
    plt.ylabel("Gene Count")
    plt.xlabel("GO Terms")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title="Region", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(Path(output_folder) / filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_grouped_tissue_analysis(region_go_mapping, output_folder):
    """
    Create grouped analysis plot (neck vs bag vs valve).
    """
    group_gene_counts = defaultdict(lambda: defaultdict(int))

    for group_name, subregions in REGION_GROUPS.items():
        for region in subregions:
            if region not in region_go_mapping:
                continue
            
            for go_term, genes in region_go_mapping[region].items():
                group_gene_counts[go_term][group_name] += len(genes)

    df = pd.DataFrame(group_gene_counts).fillna(0).T
    df = df.reindex(columns=list(REGION_GROUPS.keys()))

    plt.figure(figsize=(14,8))
    df.plot(kind='bar', stacked=True, ax=plt.gca())
    plt.title("GO Term Gene Cuonts: Neck vs Bag vs Valve")
    plt.ylabel("Gene Count")
    plt.xlabel("GO Terms")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title="Tissue Group")
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "neck_vs_bag_vs_valve.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """
    Main plotting pipeline - rune after gene lists are extracted
    """
    region_go_mapping = load_go_gene_data("go_genes_list")
    if not region_go_mapping:
        print("No GO gene data found.")
        return
    create_go_term_plots(region_go_mapping)
    print("GO plotting complete")

if __name__ == "__main__":
    main() 