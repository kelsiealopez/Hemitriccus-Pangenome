setwd("/Users/kelsielopez/Downloads")

library(readr)
library(dplyr)
library(tidyr)

df <- read_tsv("HemMar_10_haps_Oct_2025_new_faa_fix_PAV.Rtab")

## remove sex chromosome genes
gff <- read_tsv("hemMar.toga.merged.gff3",
                comment = "#", col_names = FALSE)

genes_gff <- gff %>% 
  filter(X3 == "gene") %>%
  transmute(
    scaffold   = X1,
    attributes = X9,
    Gene       = sub(".*Name=([^;]+).*", "\\1", X9)
  )

sex_chr_labels <- c("scaffold_5", "scaffold_21")
sex_chrom_genes <- genes_gff %>%
  filter(scaffold %in% sex_chr_labels) %>%
  pull(Gene) %>%
  strsplit(",") %>%
  unlist() %>%
  unique()

df <- df %>% filter(!Gene %in% sex_chrom_genes)

## haplotype columns and individuals
hap_cols   <- grep("HMRG_.*#[12]$", names(df), value = TRUE)
individuals <- unique(sub("#[12]$", "", hap_cols))



# compute raw CNV

df <- df %>%
  mutate(
    n_unique_states = apply(select(., all_of(hap_cols)), 1, \(x) n_distinct(x)),
    is_cnv = n_unique_states > 1,
    is_pav = apply(select(., all_of(hap_cols)), 1, \(x) any(x == 0))
  )

cnv_raw <- df %>% filter(is_cnv)
cat("1) Raw CNV genes (autosomes only):", nrow(cnv_raw), "\n")



# remove 0/2 CNV calls

artifact_matrix <- sapply(individuals, function(id) {
  col1 <- paste0(id, "#1")
  col2 <- paste0(id, "#2")
  (cnv_raw[[col1]] == 0 & cnv_raw[[col2]] == 2) |
    (cnv_raw[[col1]] == 2 & cnv_raw[[col2]] == 0)
})

cnv_raw$bad_any <- apply(artifact_matrix, 1, any)

cnv_no_art <- cnv_raw %>% filter(!bad_any)
cat("2) CNV genes after removing 0/2 artifacts:", nrow(cnv_no_art), "\n")



# keep cnvs present in >=2 haplotypes

cnv_in_2haps <- cnv_no_art %>%
  filter(apply(select(., all_of(hap_cols)), 1, function(x) sum(x > 0) >= 2))

cat("3) CNV genes (after artifacts) present in ≥2 haplotypes:",
    nrow(cnv_in_2haps), "\n")



# get full PAV 0/0 
# Restrict to CNVs that have at least one 0 somewhere (PAV-like)
cnv_in_2haps_pav <- cnv_in_2haps %>% filter(is_pav)

# For each individual, mark genes with both haplotypes = 0 (FULL PAV in that individual)
absent_matrix <- sapply(individuals, function(id) {
  cols <- paste0(id, c("#1", "#2"))
  apply(cnv_in_2haps_pav[, cols], 1, function(x) all(x == 0))
})

# Long-format: Gene–Individual where gene is fully absent
absent_df <- as.data.frame(absent_matrix)
absent_df$Gene <- cnv_in_2haps_pav$Gene

absent_long <- absent_df %>%
  pivot_longer(-Gene, names_to = "Individual", values_to = "Absent") %>%
  filter(Absent) %>%
  select(Gene, Individual)

cat("4a) FULL PAV gene–individual pairs among step-3 CNVs:",
    nrow(absent_long), "\n")


# load blast validated hits 
blast_hits <- read_tsv("blast_validated_hits.tsv",
                       col_types = cols(
                         Gene = col_character(),
                         Individual = col_character()
                       ))




head(blast_hits)
length(unique(blast_hits$Individual))



# drop genes that have full PAV that have a blast hit 1e-10  >= 95% identity

# $11 <= 1e-10 && $3 >= 95

# PAVs that are contradicted by BLAST (0/0 but have a good hit)
pav_with_blast <- inner_join(absent_long, blast_hits,
                             by = c("Gene", "Individual"))

genes_to_drop <- unique(pav_with_blast$Gene)

cat("4b) Genes with ≥1 FULL PAV that has a BLAST hit:",
    length(genes_to_drop), "\n")

# Final CNV set: step-3 CNVs minus those genes
cnv_final <- cnv_in_2haps %>%
  filter(!Gene %in% genes_to_drop)

cat("4c) Final CNV genes after BLAST-based PAV filtering:",
    nrow(cnv_final), "\n")


# save final lists 

writeLines(cnv_final$Gene, "final_cnv_after_blast_filter.txt")
