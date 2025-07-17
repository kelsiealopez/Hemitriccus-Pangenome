##################################################################################
# Window-based feature importance analysis with SNPs, indels, and SVs
##################################################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/feature_importance")

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(patchwork)
library(scales)
library(randomForest)
library(broom)
library(factoextra)

# Constants and labeling
pretty_labels <- c(
  recomb_rate        = "Recombination Rate",
  gene_density       = "Gene Density",
  gc_content         = "GC Content",
  dna_density        = "DNA",
  line_density       = "LINE",
  ltr_density        = "LTR",
  sine_density       = "SINE",
  satellite_density  = "Satellite",
  low_complexity_density = "Low Complexity",
  other_density      = "Other",
  simple_repeat_density  = "Simple Repeat",
  unknown_density    = "Unclassified"
)
custom_colors <- c(
  "LINE"               = alpha("#E63946", 1),
  "Unclassified"       = alpha("#F4A261", 1),
  "LTR"                = alpha("#FCDE93", 1),
  "Simple Repeat"      = alpha("#A8D5BA", 1),
  "SINE"               = alpha("#087E8B", 1),
  "Other"              = alpha("#1E3F66", 1),
  "Satellite"          = alpha("#7C1D6F", 1),
  "Low Complexity"     = alpha("#DC3977", 1),
  "DNA"                = alpha("gray", 1),
  "Recombination Rate" = "#377eb8",
  "Gene Density"       = "#356644",
  "GC Content"         = "#3fa34d"
)
pretty_order <- c(
  "Recombination Rate", "Gene Density", "GC Content", "Simple Repeat",
  "DNA", "LTR", "Unclassified", "LINE", "Low Complexity", "Satellite", "Other", "SINE"
)
vars_to_scale <- names(pretty_labels)

# Read in 1Mb window data 
sv_count    <- read_tsv("SV_count_1Mb.bed",    col_names = c('chrom', 'start', 'end', 'SV_count_label','SV_count'))     %>%
  select(chrom, start, end, SV_count)
indel_count <- read_tsv("INDEL_count_1Mb.bed", col_names = c('chrom', 'start', 'end', 'INDEL_count_label','INDEL_count')) %>%
  select(chrom, start, end, INDEL_count)
snp_count   <- read_tsv("SNP_count_1Mb.bed",   col_names = c('chrom', 'start', 'end', 'SNP_count_label','SNP_count'))     %>%
  select(chrom, start, end, SNP_count)
recomb      <- read_tsv("recomb_rate_1Mb_windows.bed",     col_names = c('chrom','start','end','recomb_rate'))
gene        <- read_tsv("gene_density_1Mb.bed",            col_names = c('chrom', 'start', 'end', 'gene_density'))
gc          <- read_tsv("gc_content.1Mb.bed",              col_names = c('chrom', 'start', 'end', 'gc_content'))
read_density <- function(densfile, base) {
  read_tsv(densfile, col_names = c('chrom','start','end',paste0(base,"_label"),paste0(base,"_density")), show_col_types=FALSE) %>%
    select(chrom, start, end, !!paste0(base,"_density"))
}
DNArep         <- read_density("density_DNA_1Mb.bed",            "dna")
LINErep        <- read_density("density_LINE_1Mb.bed",           "line")
LTRrep         <- read_density("density_LTR_1Mb.bed",            "ltr")
SINErep        <- read_density("density_SINE_1Mb.bed",           "sine")
Satellite      <- read_density("density_Satellite_1Mb.bed",      "satellite")
Low_complexity <- read_density("density_Low_complexity_1Mb.bed", "low_complexity")
Other_rep      <- read_density("density_Other_1Mb.bed",          "other")
Simple_repeat  <- read_density("density_Simple_repeat_1Mb.bed",  "simple_repeat")
Unknown_rep    <- read_density("density_Unknown_1Mb.bed",        "unknown")

# Merge all features for each variant set 
merge_preds <- function(count_tbl, y_col) {
  out <- count_tbl %>%
    left_join(recomb,             by=c('chrom','start','end')) %>%
    left_join(gene,               by=c('chrom','start','end')) %>%
    left_join(gc,                 by=c('chrom','start','end')) %>%
    left_join(DNArep,             by=c('chrom','start','end')) %>%
    left_join(LINErep,            by=c('chrom','start','end')) %>%
    left_join(LTRrep,             by=c('chrom','start','end')) %>%
    left_join(SINErep,            by=c('chrom','start','end')) %>%
    left_join(Satellite,          by=c('chrom','start','end')) %>%
    left_join(Low_complexity,     by=c('chrom','start','end')) %>%
    left_join(Other_rep,          by=c('chrom','start','end')) %>%
    left_join(Simple_repeat,      by=c('chrom','start','end')) %>%
    left_join(Unknown_rep,        by=c('chrom','start','end'))
  out %>% filter(if_all(c(y_col, vars_to_scale), ~!is.na(.)))
}
merged_sv    <- merge_preds(sv_count,    "SV_count")
merged_indel <- merge_preds(indel_count, "INDEL_count")
merged_snp   <- merge_preds(snp_count,   "SNP_count")

# Standardized regression coefficients
fit_scaled <- function(response_name, data, vars_to_scale) {
  df2 <- data %>%
    mutate(across(all_of(vars_to_scale), ~as.numeric(scale(.)), .names="s_{.col}"),
           s_y = as.numeric(scale(.data[[response_name]])))
  formula_str <- paste0("s_y ~ ", paste0("s_", vars_to_scale, collapse=" + "))
  fit <- lm(as.formula(formula_str), data=df2)
  broom::tidy(fit, conf.int=TRUE)
}
coefs_sv    <- if (nrow(merged_sv))    fit_scaled("SV_count",    merged_sv,    vars_to_scale) %>% mutate(response="SV_count") else tibble()
coefs_indel <- if (nrow(merged_indel)) fit_scaled("INDEL_count", merged_indel, vars_to_scale) %>% mutate(response="INDEL_count") else tibble()
coefs_snp   <- if (nrow(merged_snp))   fit_scaled("SNP_count",   merged_snp,   vars_to_scale) %>% mutate(response="SNP_count") else tibble()

# Prepare for plotting: format and labeling 
all_coefs_s <- bind_rows(coefs_sv, coefs_indel, coefs_snp) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = gsub("s_", "", term),
    pretty_term = pretty_labels[term],
    response = factor(response, levels = c("SNP_count", "INDEL_count", "SV_count")),
    pretty_term = factor(pretty_term, levels = rev(pretty_order)) 
  )

# Random forest feature importance 
get_rf_importance <- function(df, response, predictors) {
  rf_df <- df %>% select(all_of(c(response, predictors))) %>% na.omit()
  f <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  set.seed(123)
  rf <- randomForest(f, data=rf_df, importance=TRUE, ntree=500)
  imp <- as.data.frame(importance(rf, type=1))
  imp$term <- rownames(imp)
  imp$response <- response
  imp
}
imp_sv    <- get_rf_importance(merged_sv,    "SV_count",    vars_to_scale)
imp_indel <- get_rf_importance(merged_indel, "INDEL_count", vars_to_scale)
imp_snp   <- get_rf_importance(merged_snp,   "SNP_count",   vars_to_scale)

rf_imp_all_raw <- bind_rows(imp_sv, imp_indel, imp_snp)
# Pick correct column name automatically:
importance_col <- if ("MeanDecreaseAccuracy" %in% names(rf_imp_all_raw)) {
  "MeanDecreaseAccuracy"
} else if ("%IncMSE" %in% names(rf_imp_all_raw)) {
  "%IncMSE"
} else {
  stop("Can't find importance column in randomForest::importance() output")
}
rf_imp_all <- rf_imp_all_raw %>%
  mutate(
    pretty_term = pretty_labels[term],
    response = factor(response, levels = c("SNP_count", "INDEL_count", "SV_count")),
    pretty_term = factor(pretty_term, levels = rev(pretty_order)),
    importance = .data[[importance_col]]
  ) %>%
  arrange(pretty_term, response)

# --- Combine for scatter plots ---
combined <- merged_sv %>%
  select(chrom, start, end, all_of(vars_to_scale), SV_count) %>%
  left_join(merged_indel %>% select(chrom, start, end, INDEL_count), by=c("chrom","start","end")) %>%
  left_join(merged_snp %>% select(chrom, start, end, SNP_count), by=c("chrom","start","end"))
combined_long <- combined %>%
  pivot_longer(
    cols = c(SV_count, INDEL_count, SNP_count),
    names_to = "response",
    values_to = "count"
  )
combined_long$response <- factor(combined_long$response, levels = c("SNP_count", "INDEL_count", "SV_count"))

##################################################
########.   FINAL FIGURE PANELS    ##############
##################################################

## 1. Regression Forest Plot
p_forests <- ggplot(all_coefs_s, aes(x=estimate, y=pretty_term, color=pretty_term, shape=response)) +
  geom_point(position=position_dodge(width=0.7), size=3) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0.2, position=position_dodge(width=0.7)) +
  geom_vline(xintercept=0, lty=2, color="gray50") +
  xlab("Standardized regression coefficient") +
  ylab("") +
  theme_bw(base_size = 14) +
  scale_color_manual(
    name = "Predictor",
    values = custom_colors,
    breaks = pretty_order
  ) +
  scale_shape_manual(
    name = "Variant Type",
    values = c("SNP_count" = 15, "INDEL_count" = 17, "SV_count" = 16),
    labels = c("SNP count", "INDEL count", "SV count")
  )
print(p_forests)

## 2. Random Forest Feature Importance Bar Plot
p_rfbar <- ggplot(rf_imp_all, aes(x=importance, y=pretty_term, fill=pretty_term, group=response)) +
  geom_col(position=position_dodge(width=0.7), width=0.65, aes(alpha=response)) +
  geom_text(aes(label=sprintf("%.1f", importance), x=importance),
            position=position_dodge(width=0.7),
            vjust=0.5, hjust=-0.15, size=3,
            check_overlap=TRUE) +
  xlab("Random Forest Importance (%IncMSE)") +
  ylab("") +
  theme_bw(base_size = 14) +
  scale_fill_manual(
    name = "Predictor",
    values = custom_colors,
    breaks = pretty_order
  ) +
  scale_alpha_manual(
    name = "Variant Type",
    values = c("SNP_count" = 0.9, "INDEL_count" = 0.6, "SV_count" = 0.3),
    labels = c("SNP count", "INDEL count", "SV count")
  )
print(p_rfbar)

## 3. Patchwork: Window correlation plots
vars_ordered <- setNames(names(pretty_labels), pretty_labels)[pretty_order]
no_prop <- c("Recombination Rate", "Gene Density", "GC Content")
plot_list <- list()
for (i in seq_along(vars_ordered)) {
  pr      <- vars_ordered[i]
  pr_label <- pretty_labels[pr]
  this_color <- custom_colors[pr_label]
  x_axis_label <- if (pr_label %in% no_prop) pr_label else paste(pr_label, "proportion")
  plot_list[[i]] <- ggscatter(
    combined_long,
    x = pr, y = "count", color = "black",
    add = "reg.line",
    conf.int = TRUE,
    add.params = list(linetype = 1),
    size = 0.8, alpha = 0.4
  ) +
    facet_wrap(~response, scales = "free_y", nrow = 1) +
    stat_cor(
      method = "spearman",
      label.sep = "\n",
      aes(group = response),
      label.x.npc = 1,
      label.y.npc = 0.98,
      hjust = 1,
      vjust = 1,
      size = 3,
      show.legend = FALSE
    ) +
    labs(
      title = NULL,
      y = "Variant count",
      x = x_axis_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(color = "black", face = "bold", size = 8),
      strip.background = element_rect(fill = this_color, color = this_color),
      strip.text = element_text(face = "bold", color = "black", size = 8),
      legend.position = "none"
    )
}
print(wrap_plots(plot_list, nrow = 6, ncol = 2)) # Adjust as needed

