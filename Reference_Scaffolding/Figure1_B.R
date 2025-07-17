# Load required libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

### 1. Read and prepare data ###

# Read the RepeatMasker output (.out) file
repeat_data <- read.table("hemMar_complete_sorted_JBAT.FINAL.full_mask.wpasser.out",
                          header = FALSE, fill = TRUE, sep = "", stringsAsFactors = FALSE)

# Define column names
colnames(repeat_data) <- c("score", "perc_div", "perc_del", "perc_ins", "query",
                           "begin", "end", "left_query", "strand", "repeat",
                           "class_family", "begin_rep", "end_rep", "left_rep", "ID")

# Remove first two rows if they are not data (e.g., headers/extra lines)
repeat_data <- repeat_data[-c(1, 2), ]

# Keep just scaffolds 1:34
repeat_data <- repeat_data %>%
  filter(query %in% paste0("scaffold_", 1:34))

# Convert columns to numeric as needed
repeat_data <- repeat_data %>%
  mutate(
    begin = as.numeric(begin),
    end   = as.numeric(end)
  )

# Extract scaffold name for plotting
repeat_data <- repeat_data %>%
  mutate(scaffold = str_extract(query, "scaffold_\\d+"))

# Calculate repeat length
repeat_data <- repeat_data %>%
  mutate(length = end - begin + 1)

# Extract repeat class before the slash
repeat_data <- repeat_data %>%
  mutate(class_family_group = str_extract(class_family, "^[^/]+"))

#DNA repeats to DNA category
repeat_data <- repeat_data %>%
  mutate(broad_category = case_when(
    class_family_group %in% c("LINE", "LINE.inc", "LINE?")      ~ "LINE",
    class_family_group %in% c("DNA", "DNA.inc", "DNA?")         ~ "DNA",
    class_family_group %in% c("LTR", "LTR.inc", "LTR?")         ~ "LTR",
    class_family_group == "Low_complexity"                      ~ "Low_complexity",
    class_family_group %in% c("RC", "RC?")                      ~ "RC",
    class_family_group %in% c("Retroposon", "Retroposon?")      ~ "Retroposon",
    class_family_group %in% c("SINE", "SINE?")                  ~ "SINE",
    class_family_group == "Satellite"                            ~ "Satellite",
    class_family_group == "Segmental"                            ~ "Segmental",
    class_family_group == "Simple_repeat"                        ~ "Simple_repeat",
    class_family_group %in% c("Unknown", "Unknown.inc")          ~ "Unknown",
    class_family_group == "Unspecified"                          ~ "Unspecified",
    class_family_group %in% c("rRNA", "scRNA", "snRNA", "srpRNA", "tRNA") ~ class_family_group,
    TRUE                                                        ~ "Other"
  ))

# Calculate the total repeats and lengths per scaffold and broad category
scaffold_data <- repeat_data %>%
  group_by(scaffold, broad_category) %>%
  summarise(
    total_repeats = n(),
    total_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate total repeats per scaffold
scaffold_totals <- scaffold_data %>%
  group_by(scaffold) %>%
  summarise(total_repeats = sum(total_repeats, na.rm = TRUE))

# Join totals, calculate proportion of each repeat category for each scaffold
scaffold_data <- scaffold_data %>%
  left_join(scaffold_totals, by = "scaffold") %>%
  mutate(proportion = total_repeats.x / total_repeats.y)

# Scaffold lengths table
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

# Add scaffold lengths and coverage
scaffold_data <- scaffold_data %>%
  left_join(scaffold_lengths, by = "scaffold") %>%
  mutate(proportion_of_length = total_length / length)

# Order factors for plotting
scaffold_data <- scaffold_data %>%
  mutate(scaffold = factor(scaffold, levels = paste0("scaffold_", 1:34)))

# Add numeric/scientific scaffold labels
scaffold_data <- scaffold_data %>%
  mutate(scaffold_num = as.character(gsub("scaffold_", "", scaffold)),
         scaffold_num = case_when(
           scaffold_num == "5"  ~ "Z",
           scaffold_num == "21" ~ "W",
           TRUE                 ~ scaffold_num
         ))

scaffold_order <- c("Z", "W", setdiff(as.character(1:34), c("5", "21")))
scaffold_data <- scaffold_data %>%
  mutate(scaffold_num = factor(scaffold_num, levels = scaffold_order))

# Make sure DNA is in your color and label palettes


# this goes with other ones..... idk.
custom_colors <- c(
  "LINE"             = alpha("#E63946", 1),
  "Unclassified"     = alpha("#F4A261", 1),
  "LTR"              = alpha("#FCDE93", 1),
  "Simple Repeat"    = alpha("#A8D5BA", 1),
  "SINE"             = alpha("#087E8B", 1),
  "Other"            = alpha("#1E3F66", 1),
  "Satellite"        = alpha("#7C1D6F", 1),
  "Low Complexity"   = alpha("#DC3977", 1),
  "DNA"              = alpha("gray", 1)
)


custom_colors_alpha <- alpha(custom_colors, 0.85)

custom_labels <- c(
  "Simple_repeat"   = "Simple Repeats",
  "SINE"            = "SINE",
  "LINE"            = "LINE",
  "DNA"             = "DNA",
  "Low_complexity"  = "Low Complexity",
  "LTR"             = "LTR",
  "Unknown"         = "Unknown",
  "Other"           = "Other",
  "Satellite"       = "Satellite"
)

# --- Set factor levels for broad_category for correct plot stacking order ---
scaffold_data <- scaffold_data %>%
  mutate(broad_category = factor(broad_category,
                                 levels = c("Simple_repeat", "SINE", "LINE", "DNA", "Low_complexity", "LTR", "Unknown", "Other", "Satellite")
  ))

### FINAL PLOT ###
ggplot(scaffold_data, aes(x = scaffold_num, y = proportion_of_length, fill = broad_category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.25) +
  labs(title = "",
       x = "Scaffold",
       y = "Proportion of Scaffold",
       fill = "Repeat Type") +
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "black")
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete() +
  scale_fill_manual(values = custom_colors_alpha, labels = custom_labels)
