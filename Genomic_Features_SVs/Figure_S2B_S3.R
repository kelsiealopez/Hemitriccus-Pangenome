setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/feature_importance")

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)

# === 1. Load and prep data ===
sv_count    <- read_tsv("SV_count_1Mb.bed",    col_names = c('chrom','start','end','SV_count_label','SV_count')) %>% select(chrom, start, end, SV_count)
indel_count <- read_tsv("INDEL_count_1Mb.bed", col_names = c('chrom','start','end','INDEL_count_label','INDEL_count')) %>% select(chrom, start, end, INDEL_count)
snp_count   <- read_tsv("SNP_count_1Mb.bed",   col_names = c('chrom','start','end','SNP_count_label','SNP_count')) %>% select(chrom, start, end, SNP_count)

macro <- paste0("scaffold_", 1:17)
micro <- paste0("scaffold_", 18:34)
add_chrom_class <- function(df) {
  df %>%
    mutate(
      chrom_class = case_when(
        chrom %in% macro ~ "macro",
        chrom %in% micro ~ "micro",
        TRUE ~ NA_character_
      )
    ) %>% filter(!is.na(chrom_class))
}
snp_count_box   <- add_chrom_class(snp_count)
indel_count_box <- add_chrom_class(indel_count)
sv_count_box    <- add_chrom_class(sv_count)

# === 2. BOXPLOT GRID (facet, with p values) ===
p_sv <- ggboxplot(
  sv_count_box, x = "chrom_class", y = "SV_count", fill = "chrom_class"
) + stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("macro", "micro")),
  label = "p.format",
  label.y = max(sv_count_box$SV_count, na.rm=TRUE) * 1.05
) + labs(title = "SV", x = "")
p_indel <- ggboxplot(
  indel_count_box, x = "chrom_class", y = "INDEL_count", fill = "chrom_class"
) + stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("macro", "micro")),
  label = "p.format",
  label.y = max(indel_count_box$INDEL_count, na.rm=TRUE) * 1.05
) + labs(title = "INDEL", x = "")
p_snp <- ggboxplot(
  snp_count_box, x = "chrom_class", y = "SNP_count", fill = "chrom_class"
) + stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("macro", "micro")),
  label = "p.format",
  label.y = max(snp_count_box$SNP_count, na.rm=TRUE) * 1.05
) + labs(title = "SNP", x = "")
ggarrange(p_sv, p_indel, p_snp, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")

# === 3. SV DENSITY LINES (chromosomes 1–13, RAW SV count per 1Mb) ===
sv_windows <- sv_count %>%
  mutate(scaffold_num = as.numeric(str_extract(chrom, "[0-9]+$"))) %>%
  filter(scaffold_num >= 1, scaffold_num <= 13)

scaffold_lengths <- sv_windows %>%
  group_by(chrom) %>%
  summarize(chr_length = max(end))

sv_windows <- sv_windows %>%
  left_join(scaffold_lengths, by = "chrom") %>%
  mutate(percent = 100 * start / chr_length)

ggplot(sv_windows, aes(x = percent, y = SV_count, color = as.factor(scaffold_num), group = chrom)) +
  geom_line(linewidth = 1, alpha = 0.95) +
  scale_color_viridis_d(option = "plasma", name = "Chromosome", direction = 1) +
  labs(
    title = "",
    x = "Chromosome position (%)",
    y = "SV count per 1Mb"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.5, 0.75),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9.5)
  ) +
  guides(color = guide_legend(
    nrow = 4, byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))
