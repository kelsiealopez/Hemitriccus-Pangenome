# Plotting pie chart of whole genome repeat masking for the reference 

library(ggplot2)



# Updated repeat data including DNA with specified order
repeat_data <- data.frame(
  Category = factor(c("Nonrepetitive", "Simple_repeat", "DNA", "LTR", "Unknown", 
                      "LINE", "Low_complexity", "Satellite", "Other", "SINE"),
                    levels = c("Nonrepetitive", "Simple_repeat", "DNA", "LTR", 
                               "Unknown", "LINE", "Low_complexity", "Satellite", "Other", "SINE")),
  Percentage = c(100 - 13.86, 1.00, 0.52, 3.00, 1.07, 7.16, 0.26, 0.66, 0, 0.19)  
)



# Define custom colors, adding gray for DNA
custom_colors <- c(
  "Nonrepetitive" = alpha("white", 1),
  "LINE" = alpha("#E63946",1),       # Red
  "Unknown" = alpha("#F4A261", 1),    # Orange
  "LTR" = alpha("#FCDE93", 1),        # Yellow
  "Simple_repeat" = alpha("#A8D5BA", 1),  # More saturated green
  "SINE" = alpha("#087E8B", 1),       # Darker teal green
  "Other" = alpha("#1E3F66", 1),      # Darker blue
  "Satellite" = alpha("#7C1D6F", 1),  # Purple
  "Low_complexity" = alpha("#DC3977", 1), # Pink
  "DNA" = alpha("gray", 1)            # Gray for DNA
)

# Create a pie chart with custom colors
ggplot(repeat_data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  theme_void() +
  labs(title = "") +
  scale_fill_manual(values = custom_colors)




# fix labels to be consistent with manuscript figure 3

# Define display labels
custom_labels <- c(
  "Nonrepetitive" = "Single Copy",
  "Simple_repeat" = "Simple",
  "DNA" = "DNA",
  "LTR" = "LTR",
  "Unknown" = "Unknown",
  "LINE" = "LINE",
  "Low_complexity" = "Low complexity",
  "Satellite" = "Satellite",
  "Other" = "Other",
  "SINE" = "SINE"
)

# Plot with updated legend labels
ggplot(repeat_data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  theme_void() +
  labs(title = "") +
  scale_fill_manual(values = custom_colors, labels = custom_labels)





ggplot(repeat_data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y") +
  theme_void() +
  labs(title = "", fill = "Repeat") +  # Legend title updated to "Repeat"
  scale_fill_manual(values = custom_colors, labels = custom_labels)





# trying to add repeat content percentages instead of doing manually


