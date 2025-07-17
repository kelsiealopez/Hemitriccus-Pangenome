####### set wd and load libraries #######

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/pggb_redo/merging")
library(readr)
library(ggplot2)
library(dplyr)
library(ggsci)
library(stringr)
library(RColorBrewer)
library(viridis)
library(colorspace)
library(scales)
library(tidyr)
library(patchwork)
library(grid)

####### constants (colors, labels, broad repeat category mapping) #######

custom_colors <- c(
  "Nonrepetitive" = alpha("white", 0.8),
  "LINE" = alpha("#E63946", 0.8),
  "Unknown" = alpha("#F4A261", 0.8),
  "LTR" = alpha("#FCDE93", 0.8),
  "Multiple" = alpha("#4C9B82", 0.8),
  "Simple_repeat" = alpha("#A8D5BA", 0.8),
  "SINE" = alpha("#087E8B", 0.8),
  "Other" = alpha("#1E3F66", 0.8),
  "Satellite" = alpha("#7C1D6F", 0.8),
  "Low_complexity" = alpha("#DC3977", 0.8)
)
custom_labels <- c(
  "Nonrepetitive" = "Single copy",
  "LINE" = "LINE",
  "Unknown" = "Unclassified",
  "LTR" = "LTR",
  "Multiple" = "Multiple forms",
  "Simple_repeat" = "Simple repeat",
  "SINE" = "SINE",
  "Other" = "Other",
  "Satellite" = "Satellite",
  "Low_complexity" = "Low complexity"
)

broad_category_mapping <- c(
  "No Repeat" = "Nonrepetitive", "LINE" = "LINE", "LINE.inc" = "LINE", "LINE?" = "LINE",
  "DNA" = "Other", "DNA.inc" = "Other", "DNA?" = "Other", "LTR" = "LTR", "LTR.inc" = "LTR", "LTR?" = "LTR",
  "Low_complexity" = "Low_complexity", "Other" = "Other", "RC" = "Other", "RC?" = "Other",
  "Retroposon" = "Other", "Retroposon?" = "Other", "SINE" = "SINE", "SINE?" = "SINE",
  "Satellite" = "Satellite", "Segmental" = "Other", "Simple_repeat" = "Simple_repeat", "Unknown" = "Unknown",
  "Unknown.inc" = "Unknown", "Unspecified" = "Other", "rRNA" = "Other", "scRNA" = "Other",
  "snRNA" = "Other", "srpRNA" = "Other", "tRNA" = "Other"
)
repeat_factor_levels <- c(
  "Nonrepetitive", "Simple_repeat", "Multiple", "LTR", "Unknown",
  "LINE", "Low_complexity", "Satellite", "Other", "SINE"
)

                ########## reading and preparing data  ########## 

# Main SV data
pggb <- read_delim("pggb_variation_overlaps_final_with_max_allele_length.tab", delim = "\t")

# Separate out INS and DEL, assign negative sign to deletions
get_SV_type_df <- function(df, types, sgn=1, label) {
  dplyr::filter(df, subtype %in% types) %>%
    mutate(longest_allele_length = sgn * longest_allele_length,
           type = label)
}
SV_INS <- get_SV_type_df(pggb, c("SVINS", "SVINS_Complex"),  1, "Insertion")
SV_DEL <- get_SV_type_df(pggb, c("SVDEL", "SVDEL_Complex"), -1, "Deletion")

combined_SVs <- bind_rows(
  select(SV_INS, longest_allele_length, type),
  select(SV_DEL, longest_allele_length, type)
) %>%
  filter(abs(longest_allele_length) >= 50)  # Remove small indels

                  ######### RepeatMasker Data ########

# SV lengths (query, length)
sv_lengths <- read.table("../../repeats/SV_alleles_pggb_redo/HemMar_complete_custom_lib/length_input_query_seqs_noHemMar.txt",
                         header = FALSE, col.names = c("query", "length"))

# RepeatMasker
repeat_data <- read.table(
  "../../repeats/SV_alleles_pggb_redo/HemMar_complete_custom_lib/SV_alleles_over_50_bp.fa.noHemMar.out",
  header = FALSE, skip = 3, fill = TRUE, stringsAsFactors = FALSE)
colnames(repeat_data) <- c("SW_score", "perc_div", "perc_del", "perc_ins", "query", "query_start", "query_end", "left", "strand",
                           "matching_repeat", "repeat_class/family", "repeat_start", "repeat_end", "left_repeat", "ID")
repeat_data <- dplyr::filter(repeat_data, SW_score != "*")

repeat_agg <- repeat_data %>%
  select(query, `repeat_class/family`) %>%
  group_by(query) %>%
  summarise(repeat_classes = paste(unique(`repeat_class/family`), collapse = ","), .groups="drop")

# Join with SV length & fill NA with 'No Repeat'
sv_with_repeats <- left_join(sv_lengths, repeat_agg, by = "query") %>%
  mutate(repeat_classes = ifelse(is.na(repeat_classes), "No Repeat", repeat_classes))


map_to_broad_category <- function(repeat_classes) {
  specific_classes <- strsplit(repeat_classes, ",")[[1]]
  broad <- sapply(specific_classes, function(x) {
    x_split <- strsplit(x, "/")[[1]][1]
    ifelse(x_split %in% names(broad_category_mapping), broad_category_mapping[[x_split]], "Other")
  })
  paste0(unique(broad), collapse = ",")
}

sv_with_repeats <- sv_with_repeats %>%
  mutate(
    broad_repeat_class = sapply(repeat_classes, map_to_broad_category),
    repeat_status = ifelse(str_detect(broad_repeat_class, ","), "Multiple", broad_repeat_class)
  ) %>%
  mutate(repeat_status = factor(repeat_status, levels = repeat_factor_levels))

                  ########## Histogram plotting ########## 

plot_sv_histogram <- function(data, xlim_range, color_linewidth = 0.1, show_legend = FALSE, binwidth = 5, xlab="", ylab="Frequency") {
  ggplot(data, aes(x = length, fill = repeat_status)) +
    geom_histogram(binwidth = binwidth, color = "black", alpha = 0.85, position = "stack", linewidth = color_linewidth) +
    labs(x = xlab, y = ylab, fill = "Repeat") +
    scale_fill_manual(values = custom_colors, labels = custom_labels) +
    xlim(xlim_range[1], xlim_range[2]) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.4),
      axis.ticks.length = unit(0.12, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      legend.position = if (show_legend) c(0.95, 0.95) else "none",
      legend.justification = if (show_legend) c("right", "top") else NULL,
      legend.background = if (show_legend) element_rect(fill = alpha('white', 0.5), color = NA) else element_blank()
    )
}

                    ######### pie chart function ######### 
plot_repeat_pie <- function(count_tbl) {
  ggplot(count_tbl, aes(x = "", y = proportion, fill = repeat_status)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = custom_colors, labels = custom_labels, na.translate = FALSE)
}

                    ######### pie chart summary ######### 
repeat_status_counts <- sv_with_repeats %>%
  group_by(repeat_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(desc(proportion)) %>%
  mutate(repeat_status = factor(repeat_status, levels = repeat_factor_levels))

                          ######### make plots ######### 

# Main INS/DEL bifurcated histogram
INS_DEL_plot <- ggplot(combined_SVs, aes(x = longest_allele_length, fill = type)) +
  geom_histogram(binwidth = 50, alpha = 0.9, position = "identity") +
  scale_x_continuous("SV length (50 bp to 10 kbp)", limits = c(-10000, 10000)) +
  labs(x = "Size of Structural Variants", y = "Frequency") +
  scale_fill_manual(values = c("Insertion" = "#0b090a", "Deletion" = "#bfc0c0")) +
  theme_minimal() +
  guides(fill = guide_legend(title = "SV Type")) +
  theme(
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha('white', 0.8), color = NA)
  )

# Histograms for different length ranges
smaller_subset_plot <- plot_sv_histogram(
  sv_with_repeats, c(50, 1500), 0.05, show_legend = FALSE, binwidth=5, xlab = "SV Length (50 bp to 1500 bp)"
)

larger_subset_plot <- plot_sv_histogram(
  sv_with_repeats, c(6250, 8000), 0.05, show_legend = FALSE, binwidth=5, xlab = "SV Length (6250 bp to 8000 bp)"
)

# Pie chart
pie_chart <- plot_repeat_pie(repeat_status_counts) + theme(legend.position="none")

# PIE-INS/DEL combo: overlay pie onto INS_DEL plot
pie_grob <- ggplotGrob(pie_chart)
INS_DEL_w_pie <- INS_DEL_plot +
  annotation_custom(grob = pie_grob, xmin = -9500, xmax = -4500, ymin = 3000, ymax = 7000)

                        ######### final plots ######### 

print(INS_DEL_plot)
print(smaller_subset_plot)
print(larger_subset_plot)
print(pie_chart)
#print(INS_DEL_w_pie)
