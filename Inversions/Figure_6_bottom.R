
## load Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(grid)

## 1. Read repeat data, gene, and exon data

# Repeat data (from .bed or tab file)
repeat_data <- read.table("flanked_inversions_with_repeats_scaffold_7_INV_no_prefix.bed", sep = "\t", header = FALSE)
colnames(repeat_data) <- c("InversionScaffold", "InversionStart", "InversionEnd", "InversionSize", 
                           "FlankStart", "FlankEnd", "RepeatScaffold", "RepeatStart", "RepeatEnd", "RepeatInfo")
repeat_data$RepeatType <- sub(".*,", "", repeat_data$RepeatInfo)

# Gene data (positive/negative strands)
positive_data <- read_tsv('scaffold_7_positive_range_simple.bed', col_names = c("chrom", "start", "end", "gene"))
negative_data <- read_tsv('scaffold_7_negative_range_simple.bed', col_names = c("chrom", "start", "end", "gene"))

# RNA-gene/exon mapping
rna_to_gene_mapping <- read_tsv('rna_to_gene_mapping.txt', col_names = c("RNA_ID", "Gene_Name"))
exon_data <- read_tsv('combined_exons.txt', col_names = c("scaffold", "start", "end", "feature", "rna_id")) %>%
  inner_join(rna_to_gene_mapping, by = c("rna_id" = "RNA_ID"))

## 2. clean up for plotting

default_color <- "#bfc0c0"
highlight_color <- "#2A2E82"
# Exon coloring and position (relative stacking)
exon_data <- exon_data %>%
  mutate(
    ymin = ifelse(Gene_Name %in% positive_data$gene, -2.4, -3.1),
    ymax = ifelse(Gene_Name %in% positive_data$gene, -1.9, -2.6),
    color = ifelse(Gene_Name %in% c("RAB17", "MREG", "MLPH"), highlight_color, default_color)
  )

# One center label per gene
exon_labels <- exon_data %>%
  group_by(Gene_Name) %>%
  summarize(
    exon_start = min(start),
    exon_end = max(end),
    label_position = (exon_start + exon_end) / 2,
    .groups = "drop"
  ) %>%
  filter(!grepl('chiLan', Gene_Name) & !Gene_Name %in% c('PRLHR', 'PRLHR2')) %>%
  mutate(
    label_position = case_when(
      Gene_Name == "RAMP1" ~ label_position - 2000,
      TRUE ~ label_position
    )
  )

# Gene highlight colors
positive_data <- positive_data %>%
  mutate(
    color = case_when(
      gene %in% c("RAB17", "MREG") ~ highlight_color,
      TRUE ~ default_color
    ),
    label = ifelse(grepl("LOC", gene), "", gene)
  )
negative_data <- negative_data %>%
  mutate(
    color = case_when(
      gene == "MLPH" ~ highlight_color,
      TRUE ~ default_color
    ),
    label = ifelse(grepl("LOC", gene), "", gene),
    x_adjust = case_when(
      gene == "RAMP1" ~ 0.02,
      gene == "RBM44" ~ -0.02,
      TRUE ~ 0
    )
  )

# Categorize repeats
dna_repeats   <- filter(repeat_data, RepeatType == "DNA")
line_repeats  <- filter(repeat_data, RepeatType == "LINE")
other_repeats <- filter(repeat_data, !RepeatType %in% c("DNA", "LINE"))

# For zoom/limits
overall_start <- min(repeat_data$RepeatStart)
overall_end   <- max(repeat_data$RepeatEnd)

# Inversion coordinates 
inv_start <- unique(repeat_data$InversionStart)
inv_end   <- unique(repeat_data$InversionEnd)

## 3. Main Plot Block (modify y stacking and coloring as needed)

main_plot <- ggplot() +
  # Background for the inversion region
  geom_rect(aes(xmin = inv_start / 1e6, xmax = inv_end / 1e6, ymin = -Inf, ymax = Inf),
            fill = "#F0F0F0", alpha = 0.5, inherit.aes = FALSE) +
  geom_vline(xintercept = c(inv_start, inv_end) / 1e6, linetype = "dashed", color = "black", size = 0.3) +
  
  # DNA, LINE, OTHER repeats (edit ymin/ymax as needed for better separating bands)
  geom_rect(data = dna_repeats, aes(xmin = RepeatStart / 1e6, xmax = RepeatEnd / 1e6, ymin = -3.6, ymax = -3.4),
            fill = "black", color = "black") +
  geom_rect(data = line_repeats, aes(xmin = RepeatStart / 1e6, xmax = RepeatEnd / 1e6, ymin = -4.1, ymax = -3.9),
            fill = "black", color = "black") +
  geom_rect(data = other_repeats, aes(xmin = RepeatStart / 1e6, xmax = RepeatEnd / 1e6, ymin = -4.6, ymax = -4.4),
            fill = "black", color = "black") +
  
  # Exon rectangles
  geom_rect(data = exon_data, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = ymin, ymax = ymax),
            fill = exon_data$color, color = exon_data$color) +
  
  # Gene labels
  geom_text(data = exon_labels, aes(x = label_position / 1e6, y = ifelse(Gene_Name %in% positive_data$gene, -1.85, -2.55), label = Gene_Name),
            angle = 0, vjust = 0, size = 3, color = "black") +
  
  # (Optional) Gene rectangles and labels for positive/negative strand if you want
  geom_rect(data = positive_data, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.5, ymax = 1), fill = positive_data$color, color = positive_data$color, inherit.aes = FALSE) +
  geom_text(data = positive_data, aes(x = (start + end) / 2e6, y = 1.05, label = label), size = 3, color = "black") +
  geom_rect(data = negative_data, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -1.5, ymax = -1), fill = negative_data$color, color = negative_data$color, inherit.aes = FALSE) +
  geom_text(data = negative_data, aes(x = ((start + end) / 2e6) + x_adjust, y = -0.95, label = label), size = 3, color = "black") +
  
  # Axis/labels
  scale_x_continuous(
    limits = c(overall_start / 1e6, overall_end / 1e6),
    expand = c(0, 0),
    labels = label_number(scale = 1, suffix = "", accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = c(-4.5, -4, -3.5, -2.9, -2.2, 0.75),
    labels = c("Other", "LINE", "DNA", "-", "+", "Gene Labels")
  ) +
  labs(title = "", x = "Position (Mb)", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

print(main_plot)

# Optional: To zoom in, use limits argument in scale_x_continuous:
# zoom_start <- 31921500; zoom_end <- 32134324
# + scale_x_continuous(limits = c(zoom_start / 1e6, zoom_end / 1e6), ... )

