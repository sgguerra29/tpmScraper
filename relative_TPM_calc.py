"""
Calculate relative TPM values for spermatheca-specific gene expression analysis.
Computes relative TPM (spermatheca expression / max expression across all cell types) and identifies genes with max expression within spermatheca regions. 
"""

import pandas as pd
import os
from pathlib import Path

#mapping of spermatheca regions to functional groups
SPERMATHECA_REGION_GROUPS = {
    "Spermatheca-Uterine junction": "valve", 
    "Spermatheca neck distal": "neck", 
    "Spermatheca neck proximal": "neck", 
    "Spermatheca bag distal": "bag", 
    "Spermatheca bag proximal": "bag"
}

def load_ref_expression_data(celltype_folder):
    """
    Load expression daat from all cell types to find max expression per gene.
    """

    gene_max_expression = {}

    for filename in os.listdir(celltype_folder):
        if filename.endswith(".csv"):
            filepath = Path(celltype_folder) / filename
            df = pd.read_csv(filepath)

            for _, row in df.iterrows():
                gene_name = row['gene_short_name']
                scaled_tpm = row['scaled_TPM']

                if (gene_name not in gene_max_expression or scaled_tpm > gene_max_expression[gene_name][0]):
                    gene_max_expression[gene_name] = (scaled_tpm, filename)
        
    print(f"Loaded referenec data for{len(gene_max_expression)} genes")
    return gene_max_expression

def infer_region_group_from_filename(filename):
    """
    Determine which spermatheca group a file belongs to based on filename.
    """
    filename_lower = filename.lower().replace(" ", "_")
    for region_name, group in SPERMATHECA_REGION_GROUPS.items():
        region_key = region_name.lower().replace(" ", "_")
        if region_key in filename_lower:
            return group
    return None

def calc_relative_tpm(spermatheca_folder, ref_expression_data, output_folder):
    """
    Calculate relative TPM values for spermatheca regions.
    """
    os.makedirs(output_folder, exist_ok=True)
    spermatheca_files = set(os.listdir(spermatheca_folder))

    for filename in spermatheca_files:
        if not filename.endswith(".csv"):
            continue
        filepath = Path(spermatheca_folder) / filename
        df = pd.read_csv(filepath)

        current_group = infer_region_group_from_filename(filename)

        def calc_relative_tpm_value(row):
            """
            Calculate relative TPM for a single gene
            """
            gene_name = row['gene_short_name']
            current_tpm = row['scaled_TPM']

            ref_data = ref_expression_data.get(gene_name)
            if ref_data and ref_data[0] > 0:
                return current_tpm / ref_data[0]
            return None

        def max_in_spermatheca(row):
            """
            Check is gene's max is found in the spermatheca.
            """
            gene_name = row['gene_short_name']
            ref_data = ref_expression_data.get(gene_name)

            if ref_data:
                max_source_file = ref_data[1]
                return max_source_file in spermatheca_files
            return False

        def max_in_component(row):
            """
            Check is gene's max expression is in the same spermatheca component. 
            """
            gene_name = row['gene_short_name']
            ref_data = ref_expression_data.get(gene_name)

            if ref_data:
                max_source_file = ref_data[1]
                max_source_group = infer_region_group_from_filename(max_source_file)
                return max_source_group == current_group
            return False

        df['relative_TPM'] = df.apply(calc_relative_tpm_value, axis=1)
        df['max_in_spermatheca'] = df.apply(max_in_spermatheca, axis=1)
        df['max_in_same_component'] = df.apply(max_in_component, axis=1)

        output_path = Path(output_folder) / filename
        df.to_csv(output_path, index=False)

        print(f"Processed relative TPM for: {filename}")
        print(f" Genes with max in spermatheca: {df['max_in_spermatheca'].sum()}")
        print(f" Genes with max in the same spermatheca component: {df['max_in_same_component'].sum()}")

def main():
    """
    Main processing pipeline for relative TPM calculation
    """
    celltype_folder = "source_cell_types"
    spermatheca_folder = "spermatheca_cell_types"
    output_folder = "w_relative_TPM"

    ref_data = load_ref_expression_data(celltype_folder)
    calc_relative_tpm(spermatheca_folder, ref_data, output_folder)

    print(f"\nRelative TPM calculation complete. Results saved to: {output_folder}")

if __name__ == "__main__":
    main()