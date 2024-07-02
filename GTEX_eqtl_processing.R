# Load necessary libraries
library(dplyr)
library(biomaRt)
library(stringr)
library(readr)
library(tidyr)
library(data.table)
library(org.Hs.eg.db)
library(vroom)
library(tibble)
library(TwoSampleMR)
library(MVMR)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(ieugwasr)
library(GwasDataImport)
library(ggplot2)

# Define general output directory
output_dir <- "output/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define file paths
# These are the paths to your input eQTL data files
egenes_path <- "Artery_Aorta.v8.egenes.txt"
variant_gene_pairs_path <- "Artery_Aorta.v8.signif_variant_gene_pairs.txt"

# Step 1: Read the eQTL data files into data frames
egenes_df <- read.table(egenes_path, sep = "\t", header = TRUE)
variant_gene_pairs_df <- read.table(variant_gene_pairs_path, sep = "\t", header = TRUE)

# Step 2: Clean the data
# Remove rows with non-numeric values in 'gene_chr' column
egenes_df <- egenes_df %>%
  filter(!grepl("[^0-9]", gene_chr)) %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)),
         gene_chr = as.numeric(gsub("chr", "", gene_chr)))

# Remove rows with NA values in 'variant_id'
egenes_df <- egenes_df %>%
  filter(!is.na(variant_id))

# Rename duplicated SNP IDs
egenes_df <- egenes_df %>%
  mutate(variant_id = if_else(duplicated(variant_id), paste(variant_id, row_number(), sep = "_"), variant_id))

# Step 3: Save SNP and Probe Information
# Save SNP information to .esi file
snp_info <- egenes_df %>%
  select(variant_id, chr, variant_pos, variant_pos, ref, alt) %>%
  rename(SNP = variant_id, CHR = chr, BP = variant_pos, CM = variant_pos, A1 = ref, A2 = alt)

write.table(snp_info, file = paste0(output_dir, "eqtl3.esi"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save probe information to .epi file
probe_info <- egenes_df %>%
  select(gene_id, gene_chr, gene_start, gene_end) %>%
  rename(ID = gene_id, CHR = gene_chr, BP = gene_start, CM = gene_end)

write.table(probe_info, file = paste0(output_dir, "eqtl3.epi"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 4: Calculate Z-scores
# The Z-score standardizes the effect size (slope) by its standard error (slope_se)
egenes_df <- egenes_df %>%
  mutate(Zscore = slope / slope_se)

# Save eQTL summary statistics to .besd file
eqtl_stats <- egenes_df %>%
  select(variant_id, gene_id, Zscore) %>%
  rename(SNP = variant_id, ID = gene_id, Zscore = Zscore)

write.table(eqtl_stats, file = paste0(output_dir, "eqtl3.besd"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 5: Map Gene IDs to Gene Names using Ensembl Biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Remove version numbers from gene IDs and map to gene names
egenes_df <- egenes_df %>%
  mutate(gene_id = gsub("\\..*", "", gene_id))

gene_ids <- unique(egenes_df$gene_id)
gene_names <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                    filters = "ensembl_gene_id", 
                    values = gene_ids, 
                    mart = ensembl)

# Merge gene names with egenes_df
egenes_df <- egenes_df %>%
  left_join(gene_names, by = c("gene_id" = "ensembl_gene_id"))

# Step 6: Filter for MMP genes
mmp_genes <- egenes_df %>%
  filter(grepl("MMP", hgnc_symbol))

# Step 7: Extract Chromosome, Start, and End Positions from variant_id
mmp_genes <- mmp_genes %>%
  mutate(CHR = sapply(strsplit(variant_id, "_"), "[", 1),
         CHR = gsub("chr", "", CHR),
         START = sapply(strsplit(variant_id, "_"), "[", 2),
         END = START,
         effect_allele = str_extract(variant_id, "(?<=_)[ATGC]"),
         other_allele = str_extract(variant_id, "(?<=_)[ATGC](?=_b38)"))

# Save the filtered MMP genes
write.table(mmp_genes, file = paste0(output_dir, "MMPgenes_Artery_aorta.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

# Step 8: SNP Identification using Ensembl Biomart
snpMart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Initialize an empty dataframe for storing SNP results
rsid_df <- data.frame()
error_list <- c()

# Split dataframe into chunks for processing
n_rows <- nrow(mmp_genes)
n_chunks <- ceiling(n_rows / 50)
breaks <- c(0, seq(50, by = 50, length.out = n_chunks))
df_list <- split(mmp_genes, cut(seq_len(n_rows), breaks = breaks, labels = paste0("chunk_", seq(1, n_chunks, by = 1))))

# Retrieve SNP information for each chunk
for (i in seq_along(df_list)) {
  dataframe2 <- df_list[[i]][, c("CHR", "START", "END")]
  coords1 <- apply(dataframe2, 1, function(x) paste(x[1], x[2], x[3], sep = ":", collapse = ""))
  
  temp_result <- tryCatch({
    getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
          filters = 'chromosomal_region', 
          values = coords1, 
          mart = snpMart)
  }, error = function(e) {
    error_list <- c(error_list, i)
    return(NULL)
  })
  
  if (!is.null(temp_result)) {
    rsid_df <- rbind(rsid_df, temp_result)
  }
}

# Step 9: Merge SNP information with MMP genes
mmp_genes <- merge(mmp_genes, rsid_df[, c("refsnp_id", "chr_name", "chrom_start", "chrom_end")], 
                   by.x = c("CHR", "START", "END"), by.y = c("chr_name", "chrom_start", "chrom_end"))

# Rename columns for consistency
mmp_genes <- mmp_genes %>%
  rename(snp_col = refsnp_id, beta_col = slope, se_col = slope_se, eaf_col = maf, pval_col = pval_nominal, 
         effect_allele_col = effect_allele, other_allele_col = other_allele)

# Step 10: Save the final dataframe
saveRDS(mmp_genes, paste0(output_dir, "MMP_Artery_aorta_l.rds"))

# Load and view the processed data for verification
mmp_genes <- readRDS(paste0(output_dir, "MMP_Artery_aorta_l.rds"))
View(mmp_genes)
