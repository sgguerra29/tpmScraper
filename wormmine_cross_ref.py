"""
Cross-reference WormMine actin/myosin/caclium genes with WormSeq and CenGen expression data. 
Identifies which genes from a WormMine query (genes with GO terms containing actin, myosin, or calcium) are expressed in spermatheca regions. 
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from matplotlib_venn import venn3, venn2

def load_wormmine_data(wormmine_file):
    """
    Load and process WormMine gene data.
    """
    df = pd.read_csv(wormmine_file)
    print("WormMine file columns:", df.columns.tolist())
    print("Sample data:")
    print(df.head())

    if len(df.columns) >= 3:
        df.columns = ['gene_id', 'wbgene_id', 'go_description'] + list(df.columns[3:])

    def categorize_function(go_desc):
        """
        Categorize genes based on GO term content.
        """
        go_lower = str(go_desc).lower()
        categories = []
        if 'actin' in go_lower:
            categories.append('actin')
        if 'myosin' in go_lower:
            categories.append('myosin')
        if 'calcium' in go_lower:
            categories.append('calcium')
        return ';'.join(categories) if categories else 'other'
    
    df['functional_category'] = df['go_description'].apply(categorize_function)
    return df

def load_expression_datasets():
    """
    Load WormSeq and CenGen expression data for cross-referencing.
    """
    expression_data = {}
    wormseq_files = {
        'wormseq_neck_distal' : 'spermatheca_cell_types/Spermatheca neck distal_filtered.csv',
        'wormseq_neck_proximal' : 'spermatheca_cell_types/Spermatheca neck proximal_filtered.csv', 
        'wormseq_bag_distal': 'spermatheca_cell_types/Spermatheca bag distal_filtered.csv', 
        'wormseq_bag_proximal': 'spermatheca_cell_types/Spermatheca bag proximal_filtered.csv',
        'wormseq_sput':'spermatheca_cell_types/Spermatheca-Uterine junction_filtered.csv',
        'wormseq_merged_spermatheca': 'spermatheca_cell_types/merged_wormseq_spermatheca.csv'
    }
    cengen_files = {
        'cengen_spermatheca': 'cengen/CenGen spermatheca_filtered.csv', 
        'cengen_sput': 'cengen/CenGen sp_ut_filtered.csv'
    }

    for dataset_name, file_path in wormseq_files.items():
        if Path(file_path).exists():
            df = pd.read_csv(file_path)
            if 'gene_short_name' in df.columns:
                df = df.rename(columns={'gene_short_name': 'gene'})
            if 'max_scaled_TPM' in df.columns:
                df = df.rename(columns={'max_scaled_TPM': 'scaled_TPM'})
            expression_data[dataset_name] = df[['gene', 'scaled_TPM']].copy()
            print(f"Loaded {dataset_name}: {len(df)} genes")
        else:
            print(f"Warning: {file_path} not found")
    
    for dataset_name, file_path in cengen_files.items():
        if Path(file_path).exists():
            df = pd.read_csv(file_path)
            if 'Gene name' in df.columns:
                df = df.rename(columns={'Gene name': 'gene'})
            if 'Expression level' in df.columns:
                df = df.rename(columns={'Expression level': 'scaled_TPM'})
            expression_data[dataset_name] = df[['gene', 'scaled_TPM']].copy()
            print(f"Loaded {dataset_name}: {len(df)} genes")
        else:
            print(f"Warning: {file_path} not found")
    return expression_data

def cross_ref_genes(wormmine_df, expression_data):
    """
    Cross-reference WormMine genes with expression datasets. 
    """
    wormmine_genes = set(wormmine_df['gene_id'].dropna()) | set(wormmine_df['wbgene_id'].dropna())
    print(f"WormMine genes to cross reference: {len(wormmine_genes)}")
    crossref_results = []
    for dataset_name, expr_df in expression_data.items():
        expr_genes = set(expr_df['gene'].dropna())
        matched_genes = wormmine_genes.intersection(expr_genes)
        print(f"{dataset_name}: {len(matched_genes)} matches found")

        for gene in matched_genes:
            expr_value = expr_df[expr_df['gene'] == gene]['scaled_TPM'].iloc[0]
            wormmine_info = wormmine_df[
                (wormmine_df['gene_id'] == gene) | (wormmine_df['wbgene_id'] == gene)
            ]
            for _, row in wormmine_info.iterrows():
                crossref_results.append({
                    'gene': gene, 
                    'wbgene_id': row['wbgene_id'],
                    'go_description': row['go_description'], 
                    'functional_category': row['functional_category'],
                    'dataset': dataset_name,
                    'scaled_TPM': expr_value
                })
    return pd.DataFrame(crossref_results)

def create_crossref_visualizations(crossref_df, output_folder="wormmine_analysis"):
    """
    Create visualization for the cross-referenced data. 
    """
    os.makedirs(output_folder, exist_ok=True)
    if crossref_df.empty:
        print("No cross-referenced data found for visualization.")
        return
    create_category_expression_plot(crossref_df, output_folder)
    create_dataset_coverage_plot(crossref_df, output_folder)
    create_crossref_heatmap(crossref_df, output_folder)

def create_category_expression_plot(crossref_df, output_folder):
    """
    Create box plot of expressoin levels by functional category. 
    """
    plt.figure(figsize=(10, 6))
    plot_data = crossref_df[crossref_df['functional_category'] != 'other']
    if not plot_data.empty:
        sns.boxplot(data=plot_data, x='functional_category', y='scaled_TPM')
        plt.title('Expression Levels by Functional Category\n(Actin, Myosin, Calcium genes)')
        plt.xlabel('Functional Category')
        plt.ylabel('scaled_TPM')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(Path(output_folder) / "expression_by_category.png", dpi=300, bbox_inches='tight')
        plt.close()

def create_dataset_coverage_plot(crossref_df, output_folder):
    """
    Create bar plot showing how many genes are found in each dataset.
    """
    coverage_stats = crossref_df.groupby(['dataset', 'functional_category']).size().unstack(fill_value=0)
    plt.figure(figsize=(12, 6))
    coverage_stats.plot(kind='bar', stacked=True)
    plt.title('Gene Coverage by Dataset & Functional Category')
    plt.xlabel('Dataset')
    plt.ylabel('Number of Genes')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Functional Category', bbox_to_anchor=(1.05, 1), loc='upper left') 
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "dataset_coverage.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_crossref_heatmap(crossref_df, output_folder):
    """
    Create heatmap of gene expression across datasets.
    """
    heatmap_data = crossref_df.pivot_table(
        index='gene', 
        columns='dataset', 
        values='scaled_TPM', 
        aggfunc='mean'
    )
    if heatmap_data.empty:
        print("No data available for heatmap")
        return
    log_data = np.log10(heatmap_data.fillna(0) + 1)
    plt.figure(figsize=(10, max(6, 0.3 * len(heatmap_data))))
    sns.heatmap(
        log_data, 
        cmap='viridis', 
        cbar_kws={'label': 'log10(expression + 1)'},
        linewidths=0.1
    )
    plt.title('Expression of Actin/Myosin/Calcium Genes Across Datasets')
    plt.xlabel('Dataset')
    plt.ylabel('Gene')
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "crossref_expression_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_stats(wormmine_df, crossref_df, output_folder):
    """
    Generate summary stats and save to file. 
    """
    summary_stats = []
    total_wormmine_genes = len(wormmine_df['gene_id'].dropna().unique())
    total_crossref_genes = len(crossref_df['gene'].unique())
    summary_stats.append(f"WormMine Query Results: {total_wormmine_genes} unique genes")
    summary_stats.append(f"Genes found in expression datasets: {total_crossref_genes}")
    summary_stats.append(f"Coverage: {total_crossref_genes/total_wormmine_genes*100:.1f}%")
    summary_stats.append("")

    summary_stats.append("Breakdown by Functional Category:")
    for category in wormmine_df['functional_category'].unique():
        cat_genes_wormmine = len(wormmine_df[wormmine_df['functional_category'] == category])
        cat_genes_found = len(crossref_df[crossref_df['functional_category'] == category]['gene'].unique())
        summary_stats.append(f" {category}: {cat_genes_found}/{cat_genes_wormmine} genes found")
    summary_stats.append("")

    summary_stats.append("Breakdown by Dataset:")
    for dataset in crossref_df['dataset'].unique():
        dataset_genes = len(crossref_df[crossref_df['dataset'] == dataset]['gene'].unique())
        avg_expr = crossref_df[crossref_df['dataset'] == dataset]['scaled_TPM'].mean()
        summary_stats.append(f" {dataset}: {dataset_genes} genes, avg expression: {avg_expr:.1f}")

    summary_path = Path(output_folder) / "crossref_summary.txt"
    with open(summary_path, 'w') as f:
        f.write('\n'.join(summary_stats))
    print("Summary Statistics:")
    print('\n'.join(summary_stats))

def analyze_region_specificity(crossref_df, output_folder):
    """
    Analyze which spermatheca regions have the highest expression of each gene category.
    """
    wormseq_data = crossref_df[crossref_df['dataset'].str.contains('wormseq', na=False)]
    if wormseq_data.empty:
        print("No WormSeq data found for region analysis")
        return
    region_analysis = wormseq_data.groupby(['dataset', 'functional_category'])['scaled_TPM'].agg(['mean', 'count', 'std']).round(2)
    analysis_path = Path(output_folder) / "region_specificty_analysis.csv"
    region_analysis.to_csv(analysis_path)
    plt.figure(figsize=(12, 8))
    region_names = wormseq_data.pivot_table(
        index='functional_category',
        columns='dataset',
        values='scaled_TPM',
        aggfunc='mean'
    )
    sns.heatmap(region_names, annot=True, cmap='YlOrRd', fmt='.1f')
    plt.title('Average Expression of Functional Categories by Spermatheca Reion')
    plt.xlabel('Spermatheca Region')
    plt.ylabel('Functional Category')
    plt.tight_layout()
    plt.savefig(Path(output_folder) / "region_category_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Region specificity analysis saved to: {analysis_path}")

def main(wormmine_file="wormMine_actin_myosin_calcium.csv"):
    """
    Main pipeline for WormMine cross-referencing. 
    """
    output_folder = "wormmine_analysis"
    os.makedirs(output_folder, exist_ok=True)
    print("Loading WormMine data....")
    wormmine_df = load_wormmine_data(wormmine_file)
    print("Loading expression datasets....")
    expression_data = load_expression_datasets()
    print("Cross-referencing genes....")
    crossref_df = cross_ref_genes(wormmine_df, expression_data)
    if crossref_df.empty:
        print("No matches found between WormMine genes and expression data")
        return
    crossref_path = Path(output_folder) / "wormmine_expression_crossref.csv"
    crossref_df.to_csv(crossref_path, index=False)
    print(f"Cross-referenced data saved to: {crossref_path}")
    print("Creating visualizations....")
    create_crossref_visualizations(crossref_df, output_folder)
    print("Generating summary stats...")
    generate_summary_stats(wormmine_df, crossref_df, output_folder)
    print("Analyzing region specificity...")
    analyze_region_specificity(crossref_df, output_folder)
    print(f"WormMine cross-reference analysis complete! Results in: {output_folder}")

if __name__ == "__main__":
    main("wormMine_actin_myosin_calcium.csv")
