Spermatheca Gene Expression Analysis Pipeline:

A bioinformatics pipeline for analyzing gene expression patterns across different regions of the C.elegans spermatheca using WormSeq and CenGen datasets. 

Overview:
This pipeline performs comparitive gene expression analysis between spermatheca subregions (neck, bag, valve) to identify region-specific expression patterns and functional enrichment. 

It is important to note that while both CenGen and Wormseq offer strong exprepression data for the spermatheca and spermatheca-uterine junction, they differ slightly in their dataset's nuance. Wormseq splits the Spermatheca neck and bag each into a proximal and distal subcomponent, after finding signficant expressional differences across the spermatheca regions. Alternatively, CenGen does not differentiate its spermatheca data into neck, bag, or any of their subcomponents. Both datasets use transcripts per million (TPM) to qunantify gene expression, and therefore were straightforward in their comparisons. 

Features:
Relative TPM analysis: We deemed the relative TPM of a gene to be the scaled TPM of the gene in the given dataset divided by the max scaled TPM of that genes across all cells types published by WormSeq. We were interested in genes that were not only highly expressed in the spermatheca or valve but uniquely so. 

GO Enrichment Analysis: This pipline used g:Profiler to perform Gene Ontology enrichment for genes with a scaled TPM at or above our given threshold. The threshold for scaled TPM currently stands at 400, though it can be lowered or raised depending on the user's needs. 

Visualization Suite: The pipeline generates multiple heatmaps and comparitive plots for easier digestion of expression analysis. 
    Scaled TPM heatmaps: plot genes with scaled TPM over the given threshold by dataset region. As afformentioned, CenGen offers spermatheca vs valve, while Wormseq differentiates further. 
    GO Term Plots: For WormSeq expression data, occurence of GO term were quantified and plotted per region. Regions were also compared against one another in stacked bar plots. Please note there is currently a bug in the visualization of unique go terms per region. 
    Expression Scatterplot Across Datasets: compares genes highly expressed in both the Wormseq spermatheca (bag and valve combined) and valve and the Cengen spermatheca and valve. 
    Heatmap of Common Genes: Shows genes that were over the scaled TPM threshold in both the CenGen and WormSeq datasets as well as their scaled TPM in each of the datasets for comparison. 
    WormMine Cross Reference Expression Heatmap: Shows expression (scaled_TPM) of actin/myosin/caclium related genes in spermatheca regions. 
    WormMine Expression by Coverage: Plots Gene Coverage by dataset and functional category across spermatheca regions. 
    WormMine Expression by Region: Heatmap of the average expression of functional categories (actin, myosin, and calcium) by spermatheca region. 
    WormMine Expression by Category Plot: Plots expression levels by fucntional category (calcium, actin, myosin).

Python Dependencies:
pip install pandas numpy matplotlib seaborn pathlib gpofiler-official matplotlib-venn

Data Structure:
project_folder/
|_spermatheca_cell_types/ (5 WormSeq datasets for spermatheca and valve)
|_cengen/ (2 CenGen datasets for spermatheca and valve)
|_source_cell_types/ (WormSeq datasets for all other cell types)
|_scripts/
|_ output_folders/
|_wormMine_calcium_actin_myosin.csv/ gene list from WormMine query for genes whose GO names contain 'actin', 'myosin' and/or 'calcium. 

Complete Workflow: Run the scripts in the following order for full analysis
python relative_TPM_calc.py ... calculate relative TPM values
python enrichment_analysis.py ... filter genes with TPM over threshold for enrichment analysis
python go_analysis.py ... perform GO enrichment analysis
python extract_gene_lists_go.py ... extract GO gene lists
python go_plotting.py ... create GO visualization plots
python data_processing.py ... process and filter expression data
python merge_wormseq_spermatheca.py ... merge wormseq regions
python heatmap_generator.py ... generate expression heatmaps
python comparing_wormseq_cengen.py ... compare datasets and create final visualizations
python wormmine_cross_ref.py ... cross reference datasets with WormMine query for genes associated with actin, myosin, and/or calcium. 

Output Files:
w_relative_TPM ... relative TPM calculations with spermatheca-specific flags
filtered_for_enrichment ... gene lists filtered for GO analysis
enrichment_results ... GO enrichment analysis results
go_genes_lists ... GO terms and associated genes by region
go_term_plots ... GO analysis visualization plots
scaled_TPM_heatmap ... expression matrices and heatmap images
TPM_heatmap ... relative TPM heatmaps
combined_datasets ... cross-dataset comparison results
wormmine_analysis ... cross-refernce comparison with WormMine gene list 

Adjustable Parameters:
SCALED_TPM_THRESHOLD ... minimum expression threshold
top_n ... number of top GO terms to display
REGION_GROUPS = {} ... grouping of spermatheca sub groups (used to combined bag/neck proximal and distal)

Future Directions: 
Interestingly, CenGen offers gene expression datasets for worms across multiple development periods, as well as male and hermphroditic datasets. Future work to explore expression differences across life stages would be an interesting contination of this project. 

Filtering expression plots for what WormSeq deems "housekeeping" genes could be useful in making the output data more manageable. It may, however, also exclude genes that specifically important in the spermatheca/valve that also happen to be important in other parts of C.elegans. 

Cross referencing the datasets here with WormMine datasets for genes associated with subcomponents/functionalities of interest (involved in actin/myosin and/or calcium signaling) could help further widdle down genes of interest for experimental investigation. 



