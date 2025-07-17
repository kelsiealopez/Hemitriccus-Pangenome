setwd("/Users/kelsielopez/Desktop/Hemitriccus")

# Load required libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

# 1. Read and prepare data
repeat_data <- read.table("hemMar_complete_sorted_JBAT.FINAL.full_mask.wpasser.out",
                          header = FALSE, fill = TRUE, sep = "", stringsAsFactors = FALSE)
colnames(repeat_data) <- c("score", "perc_div", "perc_del", "perc_ins", "query",
                           "begin", "end", "left_query", "strand", "repeat",
                           "class_family", "begin_rep", "end_rep", "left_rep", "ID")

repeat_data <- repeat_data[-c(1, 2), ]  # Remove potential header lines

repeat_data <- repeat_data %>%
  filter(query %in% paste0("scaffold_", 1:34)) %>%
  mutate(
    begin = as.numeric(begin),
    end   = as.numeric(end),
    scaffold = str_extract(query, "scaffold_\\d+"),
    length = end - begin + 1,
    class_family_group = str_extract(class_family, "^[^/]+"),
    broad_category = case_when(
      class_family_group %in% c("LINE", "LINE.inc", "LINE?")      ~ "LINE",
      class_family_group %in% c("DNA", "DNA.inc", "DNA?")         ~ "DNA",
      class_family_group %in% c("LTR", "LTR.inc", "LTR?")         ~ "LTR",
      class_family_group == "Low_complexity"                      ~ "Low_complexity",
      class_family_group %in% c("RC", "RC?")                      ~ "Other",
      class_family_group %in% c("Retroposon", "Retroposon?")      ~ "Other",
      class_family_group %in% c("SINE", "SINE?")                  ~ "SINE",
      class_family_group == "Satellite"                            ~ "Satellite",
      class_family_group == "Segmental"                            ~ "Other",
      class_family_group == "Simple_repeat"                        ~ "Simple_repeat",
      class_family_group %in% c("Unknown", "Unknown.inc")          ~ "Unknown",
      class_family_group == "Unspecified"                          ~ "Other",
      class_family_group %in% c("rRNA", "scRNA", "snRNA", "srpRNA", "tRNA") ~ class_family_group,
      TRUE                                                        ~ "Other"
    )
  )

# Summarise repeat content per scaffold/category
scaffold_data <- repeat_data %>%
  group_by(scaffold, broad_category) %>%
  summarise(
    total_repeats = n(),
    total_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )

# Join totals, calculate proportions
scaffold_totals <- scaffold_data %>%
  group_by(scaffold) %>%
  summarise(total_repeats = sum(total_repeats, na.rm = TRUE))

scaffold_data <- scaffold_data %>%
  left_join(scaffold_totals, by = "scaffold") %>%
  mutate(proportion = total_repeats.x / total_repeats.y)

# -------------- SCAFFOLD LENGTHS (descending order) --------------
scaffold_lengths <- data.frame(
  scaffold = paste0("scaffold_", 1:34),
  length = c(
    190267479, 153918109, 117644792, 75196354, 72676064, 62287381,
    54550543, 36626436, 30868200, 25457091, 21659815, 21312483,
    21010199, 19919371, 18790254, 16514000, 14329195, 12025504,
    11502200, 11178082, 10580436, 7877331, 6814209, 6554688,
    6284414, 5582200, 4580250, 6731593, 2590376, 2489800,
    1889280, 1852688, 1635363, 1384001
  )
)
desc_scaffold_order <- scaffold_lengths %>%
  arrange(desc(length)) %>%
  pull(scaffold)

# Add lengths, coverage, and order
scaffold_data <- scaffold_data %>%
  left_join(scaffold_lengths, by = "scaffold") %>%
  mutate(proportion_of_length = total_length / length,
         scaffold = factor(scaffold, levels = desc_scaffold_order),
         scaffold_num = as.numeric(gsub("scaffold_", "", scaffold)))

# Remove rows for categories not present for each scaffold
scaffold_data <- scaffold_data %>% filter(!is.na(broad_category))

# --- Color and label palettes ---
custom_colors <- c(
  "LINE"             = alpha("#E63946", 0.95),
  "Unknown"          = alpha("#F4A261", 0.95),
  "LTR"              = alpha("#FCDE93", 0.95),
  "Simple_repeat"    = alpha("#A8D5BA", 0.95),
  "SINE"             = alpha("#087E8B", 0.95),
  "Other"            = alpha("#1E3F66", 0.95),
  "Satellite"        = alpha("#7C1D6F", 0.95),
  "Low_complexity"   = alpha("#DC3977", 0.95),
  "DNA"              = alpha("gray", 0.95)
)
custom_colors_alpha <- alpha(custom_colors, 0.95)

custom_labels <- c(
  "Simple_repeat"   = "Simple Repeats",
  "SINE"            = "SINE",
  "LINE"            = "LINE",
  "DNA"             = "DNA",
  "Low_complexity"  = "Low Complexity",
  "LTR"             = "LTR",
  "Unknown"         = "Unclassified",
  "Other"           = "Other",
  "Satellite"       = "Satellite"
)

scaffold_data <- scaffold_data %>%
  mutate(broad_category = factor(
    broad_category,
    levels = c("Simple_repeat", "SINE", "LINE", "DNA", "Low_complexity", "LTR", "Unknown", "Other", "Satellite")
  ))

# safe to drop NAs I guess
scaffold_data <- scaffold_data %>% filter(!is.na(broad_category))



scaffold_data <- scaffold_data %>%
  mutate(broad_category = factor(
    broad_category,
    levels = c("Simple_repeat", "DNA", "LTR", "Unknown", "LINE", "Low_complexity", "Satellite", "Other", "SINE")
  ))

# --- PLOT ---
ggplot(scaffold_data, aes(x = scaffold, y = proportion_of_length, fill = broad_category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.25) +
  labs(title = "",
       x = "Scaffold",
       y = "Proportion",
       fill = "Repeat Type") +
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "black")
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = function(scaf) as.numeric(gsub("scaffold_", "", scaf))) +
  scale_fill_manual(values = custom_colors_alpha, labels = custom_labels)





