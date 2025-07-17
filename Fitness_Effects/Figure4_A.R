# looking at ratios of SV/SNP and INDEL/SNP of Hemitriccus margaritaceiventer (HM) , House Finch (HF) , and
# Scrub Jays (Aphelocoma caerulescens (AC); A. woodhouseii (AW); A. insularis (AI) from Fang & Edwards 2024, and Edwards et al. 2025

library(ggplot2)


# Create a data frame to hold your data
data <- data.frame(
  Comparison = rep(c("SV/SNP", "INDEL/SNP"), each = 5),
  Type = rep(c("HF", "AW", "AC", "AI", "HM"), times = 2),
  Ratio = c(0.031, 0.019, 0.022, 0.099, 0.015, 0.19, 0.163, 0.186, 0.451, 0.258)
)

# Reorder the Type factor according to the desired order
data$Type <- factor(data$Type, levels = c("HM", "HF", "AW", "AC", "AI"))

# Ensure the Comparison factor places SV/SNP first, then INDEL/SNP
data$Comparison <- factor(data$Comparison, levels = c("SV/SNP", "INDEL/SNP"))

# Create the plot with reduced bin width and black outlines
plot <- ggplot(data, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, 
           color = "black", size = 0.2) +  # Add black outline with thickness 0.3
  labs(title = "", 
       x = "", 
       y = "Ratio", 
       fill = "") +
  theme_classic() +
  theme(
    legend.position = c(0.05, 1),  # Position at top left inside the plot
    legend.justification = c("left", "top"),  # Anchor legend position
    legend.background = element_rect(fill = "white", color = NA)  # Optional: background for readability
  ) +
  scale_fill_manual(values = c("HM" = "#81B29A", "HF" = "#B1736C", "AW" = "#9BADE9", "AC" = "#D1D1D1", "AI" = "#385FE9"))

# Display the plot
print(plot)
