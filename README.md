
# GTEx_MMP_Analysis

This repository contains scripts and data for analyzing matrix metalloproteinase (MMP) genes using eQTL data from the GTEx (Genotype-Tissue Expression) Project.

## Overview

This project aims to identify significant eQTLs for MMP genes across different tissues using the GTEx Analysis V8 release. The analysis includes data cleaning, mapping gene IDs to gene names, calculating Z-scores for eQTLs, and retrieving SNP information.

## Data Source

The data used in this project is obtained from the GTEx Portal. Specifically, we use the single-tissue cis-eQTL data:

- `*.egenes.txt.gz`: Contains data for all genes tested; to obtain the list of eGenes, select rows with 'qval' ≤ 0.05.
- `*.signif_variant_gene_pairs.txt.gz`: Contains significant variant-gene associations based on permutations.

## Prerequisites

The following R packages are required to run the scripts in this repository:

- dplyr
- biomaRt
- stringr
- readr
- tidyr
- data.table

You can install these packages using the following command:

```r
install.packages(c("dplyr", "biomaRt", "stringr", "readr", "tidyr", "data.table"))
```

## Directory Structure

- `scripts/`: Contains the main analysis script.
- `data/`: Directory to store input data files.
- `output/`: Directory to store output files.

## Usage

1. **Download the Data**: Download the required data files from the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl).

2. **Place the Data**: Place the downloaded data files (`*.egenes.txt.gz` and `*.signif_variant_gene_pairs.txt.gz`) into the `data/` directory.

3. **Run the Script**: Execute the R script located in the `scripts/` directory to perform the analysis.

## Step-by-Step Guide

1. **Load Necessary Libraries**: The script begins by loading the required R packages.
2. **Define File Paths**: Paths to input data files and output directory are specified.
3. **Read and Clean Data**: The script reads the eQTL data files, cleans the data by removing non-numeric chromosome values and handling duplicated SNP IDs.
4. **Save SNP and Probe Information**: SNP and probe information is extracted and saved in .esi and .epi files, respectively.
5. **Calculate Z-scores**: Z-scores for the eQTLs are calculated to standardize the effect sizes.
6. **Map Gene IDs to Gene Names**: Gene IDs are mapped to gene names using the Ensembl BioMart database.
7. **Filter for MMP Genes**: The script filters the data to include only MMP genes.
8. **Extract Chromosome, Start, and End Positions**: These positions are extracted from the variant IDs for further analysis.
9. **SNP Identification**: SNP information is retrieved from the Ensembl BioMart database.
10. **Merge and Save Data**: The SNP information is merged with the filtered MMP genes data and saved as an RDS file.

## Example

```r
# Load necessary libraries
library(dplyr)
library(biomaRt)
library(stringr)
library(readr)
library(tidyr)
library(data.table)

# Define general output directory
output_dir <- "output/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define file paths
egenes_path <- "data/Artery_Aorta.v8.egenes.txt"
variant_gene_pairs_path <- "data/Artery_Aorta.v8.signif_variant_gene_pairs.txt"

# Step-by-step script as provided in the main analysis script
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

If you have any questions or need further assistance, please contact:

- Simranjit Kaur Kang
- simivk1991@gmail.com

```

This README file provides an overview of the project, instructions on how to set up the environment, and a step-by-step guide on how to run the analysis script. Feel free to customize it further based on your specific requirements.
