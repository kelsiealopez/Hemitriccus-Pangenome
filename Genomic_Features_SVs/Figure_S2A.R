# Plot Recomb Rate in 1Mb winodws across chromosomes

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scico)    # for color scale

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing")

# Load data (replace with your path if needed)
data <- read.table("cleaned_all_chroms_combined_no_missing.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preprocessing: focus on macrochromosomes
window_size <- 1e6   # 1 Mb
chrom.max   <- 13    # or set to 17 for all macrochromosomes

# Assign scaffold/chromosome number
data <- data %>%
  mutate(
    scaffold_num = as.numeric(str_extract(chrom, "[0-9]+$"))
  ) %>%
  filter(scaffold_num >= 1, scaffold_num <= chrom.max)

# Get chromosome lengths
scaffold_lengths <- data %>%
  group_by(chrom) %>%
  summarise(chr_length = max(end, na.rm = TRUE), .groups = "drop")

# Bin records to 1Mb windows
data_binned <- data %>%
  mutate(bin = floor(start / window_size) * window_size)

# Aggregate in each chrom/bin (mean recomb rates and bin start/end)
agg_data <- data_binned %>%
  group_by(chrom, scaffold_num, bin) %>%
  summarise(
    avg_recombRate = mean(recombRate, na.rm = TRUE),
    avg_CI95LO     = mean(CI95LO, na.rm = TRUE),
    avg_CI95HI     = mean(CI95HI, na.rm = TRUE),
    win_start      = min(start, na.rm = TRUE),
    win_end        = max(end, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(scaffold_lengths, by = "chrom") %>%
  mutate(percent = 100 * (win_start + win_end) / 2 / chr_length)

#  Main line plot (recombination rate along each chrom, colored by chrom 
ggplot(agg_data, aes(x = percent, y = avg_recombRate, color = as.factor(scaffold_num), group = chrom)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  scale_color_viridis_d(option = "plasma", name = "Chromosome", direction = 1) +
  labs(
    title = "",
    x = "Chromosome position (%)",
    y = "Average recombination rate per 1Mb"
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
    nrow = 4,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))

# Panel: physical position, recomb+CI, across all chroms
# To get classic linear plots of each chrom with recomb + CI, faceted panels
ggplot(agg_data, aes(x = (win_start + win_end)/2, y = avg_recombRate)) +
  geom_line(color = "blue", linewidth = 0.9) +
  geom_ribbon(aes(ymin = avg_CI95LO, ymax = avg_CI95HI), fill = "lightblue", alpha = 0.5) +
  facet_wrap(~ chrom, scales = "free_x") +
  labs(
    title = "Recombination rate (1Mb bins, by chromosome)",
    x = "Genomic position (bp)",
    y = "Avg. recombination rate"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"))



########################################################################################################




# Load dadta
data <- read.table("cleaned_all_chroms_combined_no_missing.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

window_size <- 1e6

# Assign scaffold numbers, filter for 1-13, and get per-chromosome lengths
data <- data %>%
  mutate(scaffold_num = as.numeric(str_extract(chrom, "[0-9]+$"))) %>%
  filter(scaffold_num >= 1, scaffold_num <= 13)

scaffold_lengths <- data %>%
  group_by(chrom) %>%
  summarise(chr_length = max(end, na.rm = TRUE), .groups = "drop")

# Bin to 1Mb windows and aggregate
agg_data <- data %>%
  mutate(bin = floor(start / window_size) * window_size) %>%
  group_by(chrom, scaffold_num, bin) %>%
  summarise(
    avg_recombRate = mean(recombRate, na.rm = TRUE),
    win_start = min(start, na.rm = TRUE),
    win_end = max(end, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(scaffold_lengths, by = "chrom") %>%
  mutate(percent = 100 * (win_start + win_end)/2 / chr_length)

# Final plot - matches your favorite:
ggplot(agg_data, aes(x = percent, y = avg_recombRate, color = as.factor(scaffold_num), group = chrom)) +
  geom_line(linewidth=1, alpha=0.9) +
  scale_color_viridis_d(option = "plasma", name = "Chromosome (> 20Mb)", direction = 1) +
  labs(
    title = "",
    x = "Chromosome position (%)",
    y = "Average recombination rate per 1Mb"
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
    nrow = 4,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))
