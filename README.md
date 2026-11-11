# River biofilm microeukaryotic community analysis

R scripts to process and analyse microeukaryotic 18S rRNA gene sequencing data from river biofilms.  

### Initial processing
Raw reads have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number XXXXXXXXX
- Primers are removed from raw sequencing reads with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) using DADA2/cutadapt_script.R
- Reads are processed according to the [DADA2 workflow for big data](https://benjjneb.github.io/dada2/bigdata.html) using DADA2/dada2_script.R
- Taxonomy is assigned according to the [PR2 v5.1.1](https://pr2-database.org/) reference database using DADA2/merge_and_taxonomy_script.R
- Rarefaction and filtering is performed using DADA2/filter_and_rarefy_script.R

### Analysis and visualisation
Outputs (available in the [Zenodo archive](https://doi.org/10/5281/zenodo.24136801)) can then be analysed with the provided scripts to:
   * Assign trophic classes using [EukFunc](https://github.com/lentendu/EukFunc) and analysis_scripts/taxonomy_and_EukFunc_script.R
   * Assess the diversity, community structure, and composition of the core vs non-core communities using analysis_scripts/community_structure_script.R
   * Run ordination and correlation analyses to identify environmental drivers using analysis_scripts/ordination_script.R
   * Construct co-occurrence networks using analysis_scripts/seasonal_network_script.R
   * Identify environmental thresholds and indicator taxa using analysis_scripts/titan_run.R, visualised using analysis_scripts/titan_plots_script.R

### Citation
LINK
