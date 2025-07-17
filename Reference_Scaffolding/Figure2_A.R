#  HemMar Karyotype Rideogram  Script for figure2

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/pggb_redo/merging")
library(RIdeogram)
library(dplyr)
library(readr)
library(svglite)
library(rsvg)
library(grid)

#  1. Load and format karyotype and SV density data 

hemMar_karyotype <- read.csv("chromosome_data.csv")[, 1:3]
colnames(hemMar_karyotype) <- c("Chr", "Start", "End")
hemMar_karyotype$Chr <- as.numeric(gsub("HemMar#1#scaffold_", "", hemMar_karyotype$Chr))

# Load and format SV density as value per window
SV_density <- read.csv("variant_densities.csv")
colnames(SV_density) <- c("Chr", "Start", "End", "VariantCount", "Density", "ChromNum")
SV_density <- SV_density[, c("Chr", "Start", "End", "Density")]
colnames(SV_density) <- c("Chr", "Start", "End", "Value")
SV_density$Chr <- as.numeric(SV_density$Chr)

# 1a. (Optional) Cap "variant_count" at 250 for visualization 
if ("variant_count" %in% colnames(SV_density)) {
  SV_density$Value <- pmin(SV_density$Value, 250)
}

#  1b. Ensure SV_density doesn't extend past chr ends
SV_density <- SV_density %>%
  left_join(hemMar_karyotype, by = "Chr", suffix = c(".density", ".chromosome")) %>%
  mutate(End = pmin(End.density, End.chromosome)) %>%
  select(Chr, Start = Start.density, End, Value)

# 2. Make the Rideogram with SV density as a gray heatmap 
colorset_gray <- c("#ffffff", "#343434")
ideogram(karyotype = hemMar_karyotype, overlaid = SV_density, colorset1 = colorset_gray)
convertSVG("chromosome.svg", device = "png")

# 3. Optional: add gene density heatmap overlay (requires GFFex)
# Write karyotype table for GFFex
write.table(hemMar_karyotype, "hemMar_karyotype.txt", sep = "\t", row.names = FALSE, quote = FALSE)
gene_density <- GFFex(
  input = "/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar_without_scaffold_prefix.gff3",
  karyotype = "hemMar_karyotype.txt",
  feature = "gene", window = 1e6
)

# Plot gene density as a *separate* overlaid heatmap or as a line
# Option 1: gene density as second heatmap
ideogram(
  karyotype = hemMar_karyotype,
  overlaid = gene_density,
  label = SV_density,
  label_type = "heatmap",
  colorset2 = colorset_gray
)
convertSVG("chromosome.svg", device = "png")

# Option 2: gene density as a line overlay (orange line)
gene_density$Color <- "fc8d62"  # orange hex
ideogram(
  karyotype = hemMar_karyotype,
  overlaid = SV_density,
  label = gene_density,
  label_type = "line",
  colorset1 = colorset_gray
)
convertSVG("chromosome.svg", device = "png")

# 4. View the SVG (chromosome.svg) directly as a raster in R
if (file.exists("chromosome.svg")) {
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate/locate SVG file.")
}
