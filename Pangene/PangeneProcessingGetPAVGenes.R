
setwd("/Users/kelsielopez/Downloads")

library(readr)
library(dplyr)

#Load data

df <- read_tsv("/Users/kelsielopez/Downloads/HemMar_10_haps_Oct_2025_new_faa_fix_PAV.Rtab")

# remove sex chromosome genes 

gff <- read_tsv("/Users/kelsielopez/Downloads/hemMar.toga.merged.gff3",
                comment = "#", col_names = FALSE)

genes_gff <- gff %>% filter(X3 == "gene") %>%
  transmute(
    scaffold = X1,
    attributes = X9,
    Gene = sub(".*Name=([^;]+).*", "\\1", X9)
  )

sex_chr_labels <- c("scaffold_5", "scaffold_21")
sex_chrom_genes <- genes_gff %>%
  filter(scaffold %in% sex_chr_labels) %>%
  pull(Gene) %>%
  strsplit(",") %>%
  unlist() %>%
  unique()

# remove sex chrom gene rows
df <- df %>% filter(!Gene %in% sex_chrom_genes)

# CNV/PAV Calls and Filtering 

hap_cols <- grep("HMRG_.*#[12]$", names(df), value = TRUE)
individuals <- unique(sub("#[12]$", "", hap_cols))

# CNV/PAV flag per gene
df <- df %>%
  mutate(
    n_unique_states = apply(select(., all_of(hap_cols)), 1, \(x) n_distinct(x)),
    is_cnv = n_unique_states > 1,
    is_pav = apply(select(., all_of(hap_cols)), 1, \(x) any(x == 0))
  )

# remove genes with 0/2 (or 2/0) artifact in any individual as these are likely due to misassemblies
artifact_matrix <- sapply(individuals, function(id) {
  col1 <- paste0(id, "#1")
  col2 <- paste0(id, "#2")
  (df[[col1]] == 0 & df[[col2]] == 2) | (df[[col1]] == 2 & df[[col2]] == 0)
})
df$bad_any <- apply(artifact_matrix, 1, any)

df_cnv_artifact_filtered <- df %>% filter(is_cnv, !bad_any)

#
# editing 12/3/2025
#



# --- Genes entirely absent (both haplotypes = 0) in at least one individual,
#     at the stage df_cnv_artifact_filtered (CNVs after 2/0 artifact filter) ---

# Restrict to CNVs with PAV (has at least one 0 somewhere)
df_cnv_artifact_pav <- df_cnv_artifact_filtered %>% filter(is_pav)

# For each individual, mark genes where both haplotypes are 0 (entirely absent)
absent_matrix_artifact <- sapply(individuals, function(id) {
  cols <- paste0(id, c("#1", "#2"))
  apply(df_cnv_artifact_pav[, cols], 1, function(x) all(x == 0))
})

# Total number of genes entirely absent in at least one individual
genes_absent_any_artifact <- df_cnv_artifact_pav$Gene[apply(absent_matrix_artifact, 1, any)]
cat("CNV+PAV genes (after 2/0 artifact filter) entirely absent in at least one individual:",
    length(genes_absent_any_artifact), "\n")

# Which individuals they are absent in
#  per-individual lists
absent_by_individual <- lapply(seq_along(individuals), function(i) {
  df_cnv_artifact_pav$Gene[absent_matrix_artifact[, i]]
})
names(absent_by_individual) <- individuals

# print summary
for (id in individuals) {
  cat(sprintf("Genes entirely absent in %s: %d\n",
              id, length(absent_by_individual[[id]])))
}

# long-format table where gene is absent in an individual
library(tidyr)

absent_df_artifact <- as.data.frame(absent_matrix_artifact)
absent_df_artifact$Gene <- df_cnv_artifact_pav$Gene

absent_long_artifact <- absent_df_artifact %>%
  pivot_longer(-Gene, names_to = "Individual", values_to = "Absent") %>%
  filter(Absent) %>%
  select(Gene, Individual)

# Save table 
write.table(absent_long_artifact,
            file = "cnv_artifact_filtered_genes_absent_by_individual.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# 




