#  HemMar Karyotype Rideogram  Script for figure2
# March 31st, 2025
# Set working directory where your PAF file is located
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/pggb_redo/merging")


# trying to improve karyotype plot with Rideogram 

library(RIdeogram)
library(pafr)
library(dplyr)
library(svglite)
library(rsvg)

# using data from ~/pggb_redo/merging/SV_density_karyotype_good_250331.R

# Load chromosome_data
chromosome_data <- read.csv("chromosome_data.csv")

# Load variant_densities
variant_densities <- read.csv("variant_densities.csv")
head(variant_densities)

min(variant_densities$density)
mean(variant_densities$density)
max(variant_densities$density)

# mean 9.950763e-05

# https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html


####### converting karyotype data to be more similar to Rideogram ###
# Rename the data frame to follow this tutorial 
hemMar_karyotype <- chromosome_data

# Select only the first three columns
hemMar_karyotype <- hemMar_karyotype[, 1:3]  # Alternatively, hemMar_karyotype <- hemMar_karyotype[, c("scaffold", "start", "end")]

# Rename the columns
colnames(hemMar_karyotype) <- c("Chr", "Start", "End")

# Remove the 'HemMar#1#scaffold_' prefix
hemMar_karyotype$Chr <- gsub("HemMar#1#scaffold_", "", hemMar_karyotype$Chr)

# Display the resulting data frame
print(hemMar_karyotype)


###### converting varient data to be more similar to Rideogram #####

# Create a copy of the data frame and name it SV_density
SV_density <- variant_densities

# Rename the columns
colnames(SV_density) <- c("Chr", "Start", "End", "VariantCount", "Density", "ChromNum")

# Convert the 'Start' and 'End' columns back to regular numeric format
SV_density$Start <- as.numeric(SV_density$Start)
SV_density$End <- as.numeric(SV_density$End)

# Select and rename the columns 'Chr', 'Start', 'End', and 'Density' to 'Chr', 'Start', 'End', and 'Value'
SV_density <- SV_density[, c("Chr", "Start", "End", "Density")]
colnames(SV_density) <- c("Chr", "Start", "End", "Value")

# Display the head of the transformed data frame
head(SV_density)



#ideogram(karyotype = hemMar_karyotype)
#convertSVG("chromosome.svg", device = "png")


#ideogram(karyotype = hemMar_karyotype, overlaid = SV_density)
#convertSVG("chromosome.svg", device = "png")

# fix the color 
ideogram(karyotype = hemMar_karyotype, overlaid = SV_density, colorset1 = c("#f8f9fa", "#e9ecef", "#dee2e6", "#ced4da", "#adb5bd", "#6c757d" ,"#495057", "#343a40", "#212529"))
convertSVG("chromosome.svg", device = "png")


# Convert and visualize the SVG as needed
# Ensure the file path is correct before attempting conversion
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # To view the SVG plot as a raster image directly:
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}


min(variant_densities$density)
mean(variant_densities$density)
max(variant_densities$density)



#######try to modify so anything over 250 is dark 

# Cap the variant_count at 250
variant_densities$variant_count <- ifelse(variant_densities$variant_count > 250, 250, variant_densities$variant_count)

# Recalculate the density using the capped variant_count
# Assuming the density was originally calculated as variant_count divided by the window size - which might be 1,000,000 based on the range:
variant_densities$density <- variant_densities$variant_count / (variant_densities$next_window - variant_densities$window)

# Now transform the data frame as before to prepare for the Rideogram
SV_density <- variant_densities

# Rename columns to match the Rideogram format
colnames(SV_density) <- c("Chr", "Start", "End", "VariantCount", "Density", "ChromNum")

# Convert the 'Start' and 'End' columns back to regular numeric format
SV_density$Start <- as.numeric(SV_density$Start)
SV_density$End <- as.numeric(SV_density$End)

# Select and rename the columns 'Chr', 'Start', 'End', and 'Density' to 'Chr', 'Start', 'End', and 'Value'
SV_density <- SV_density[, c("Chr", "Start", "End", "Density")]
colnames(SV_density) <- c("Chr", "Start", "End", "Value")

# Display the head of the transformed data frame
head(SV_density)


# fix the color 
ideogram(karyotype = hemMar_karyotype, overlaid = SV_density, colorset1 = c("#f8f9fa", "#e9ecef", "#dee2e6", "#ced4da", "#adb5bd", "#6c757d" ,"#495057", "#343a40", "#212529"))
convertSVG("chromosome.svg", device = "png")

?ideogram()

# Convert and visualize the SVG as needed
# Ensure the file path is correct before attempting conversion
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # To view the SVG plot as a raster image directly:
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}



# Define 18 colors to create a gradient from white to very dark gray
colorset1 <- c(
  "#ffffff", "#f8f8f8", "#f1f1f1", "#eaeaea", "#e3e3e3",
  "#dcdcdc", "#d5d5d5", "#cecece", "#c7c7c7", "#c0c0c0",
  "#b9b9b9", "#b2b2b2", "#ababab", "#a4a4a4", "#9d9d9d",
  "#969696", "#8f8f8f", "#888888", "#818181", "#7a7a7a",
  "#737373", "#6c6c6c", "#656565", "#5e5e5e", "#575757",
  "#505050", "#494949", "#424242", "#3b3b3b", "#343434"
)

# Use the color palette in your ideogram
ideogram(
  karyotype = hemMar_karyotype, 
  overlaid = SV_density, 
  colorset1 = colorset1
)

# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}




# Define 2 colors to create a gradient from white to very dark gray
colorset1 <- c(
  "#ffffff", "#343434"
)

# Use the color palette in your ideogram
ideogram(
  karyotype = hemMar_karyotype, 
  overlaid = SV_density, 
  colorset1 = colorset1
)

# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}


hemMar_karyotype
SV_density


####i need to fix the end of the SV_density data because it is extending over the karyotype data

library(dplyr)

# Ensure that both 'Chr' columns are of the same type
SV_density$Chr <- as.numeric(as.character(SV_density$Chr))
hemMar_karyotype$Chr <- as.numeric(as.character(hemMar_karyotype$Chr))

# Join `SV_density` with `hemMar_karyotype` to get the chromosome lengths
SV_density_adjusted <- SV_density %>%
  left_join(hemMar_karyotype, by = "Chr", suffix = c(".density", ".chromosome"))

# Adjust the `End.density` values not to exceed `End.chromosome`
SV_density_adjusted <- SV_density_adjusted %>%
  mutate(End = ifelse(End.density > End.chromosome, End.chromosome, End.density)) %>%
  select(Chr, Start = Start.density, End, Value)

# Display the head of the adjusted data frame
head(SV_density_adjusted)


# Define 2 colors to create a gradient from white to very dark gray
colorset1 <- c(
  "#ffffff", "#343434"
)

# Use the color palette in your ideogram
ideogram(
  karyotype = hemMar_karyotype, 
  overlaid = SV_density_adjusted, 
  colorset1 = colorset1
)


# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}



head(hemMar_karyotype)




# Use the color palette in your ideogram
ideogram(
  karyotype = hemMar_karyotype, 
  overlaid = SV_density_adjusted, 
  colorset1 = colorset1,
  width=300
)

# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")
(200000000*516) / 190267479


(190000000*516) / 190267479



(150000000*516 ) / 190267479


(100000000*516 ) / 190267479

(175000000*516 ) / 190267479

(145000000*516 ) / 190267479

(125000000*516 ) / 190267479


(50000000*516 ) / 190267479


(25000000*516 ) / 190267479


(75000000*516 ) / 190267479




# saves to ~/pggb_redo/merging

# Write the data frame to a text file
write.table(hemMar_karyotype, "hemMar_karyotype.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# trying to plot gene density

gene_density <- GFFex(input = "/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar_without_scaffold_prefix.gff3", karyotype = "hemMar_karyotype.txt", feature = "gene", window = 1000000)

sum(gene_density$Value)



# Use the color palette in your ideogram
ideogram(
  karyotype = hemMar_karyotype, 
  overlaid = gene_density
)

# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}





ideogram(karyotype = hemMar_karyotype, overlaid = gene_density, label = SV_density_adjusted, label_type = "heatmap", colorset2 = c("#f7f7f7", "#343434")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.


# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}


head(gene_density)

### trying to make it a line instead 

# Assuming gene_density is already defined
# Add a new column 'Color' with the specified color code
gene_density_with_color <- gene_density
gene_density_with_color$Color <- "fc8d62"

# Display the head of the new data frame
head(gene_density_with_color)



ideogram(karyotype = hemMar_karyotype, overlaid = SV_density_adjusted, label = gene_density_with_color, label_type = "line", colorset1 = c("#ffffff", "#343434"))

# Convert the SVG file to a PNG
convertSVG("chromosome.svg", device = "png")

# Convert and visualize the SVG as needed
if (file.exists("chromosome.svg")) {
  convertSVG("chromosome.svg", file = "chromosome", device = "png")
  
  # View the SVG plot as a raster image directly
  svg_content <- rsvg::rsvg("chromosome.svg")
  grid::grid.raster(svg_content)
} else {
  message("Error: unable to generate SVG file, check data inputs.")
}

