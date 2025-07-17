# Set working directory
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/inversions_redo")

# Libraries
library(RIdeogram)
library(dplyr)
library(rsvg)
library(grid)

#  1. Define karyotype and synteny (basic version, then post-shift) 

# Manual karyotype data: two chromosomes or species
karyotype_manual <- data.frame(
  Chr     = c(7, 7),
  Start   = c(1, 1),
  End     = c(1313376, 1317762),
  fill    = c("FFFFFF", "FFFFFF"),
  species = c("hemMar", "6433"),
  size    = c(12, 12),
  color   = c("FFFFFF", "FFFFFF")
)

# Manual synteny blocks: with a known inversion in block 2 (we'll flip Start_2, End_2 on row 2)
synteny_manual <- data.frame(
  Species_1 = rep(1, 3),
  Start_1   = c(1, 742934, 1049212),
  End_1     = c(750018, 1048957, 1313376),
  Species_2 = rep(1, 3),
  Start_2   = c(1, 748486, 1048153),      # Note: row 2 needs (End_2, Start_2) swap for inversion, later
  End_2     = c(753255, 1050206, 1317762),
  fill      = c("DCDDDE", "81b29a", "DCDDDE") # pale/green/pale, for contrast
)

#  2.  Shift coordinates as needed for alignment 
# To "move down" by 400,000 bases, subtract (ensure min 1)
karyotype_manual <- karyotype_manual %>%
  mutate(
    Start = pmax(Start - 400000, 1),
    End = pmax(End - 400000, 1)
  )
synteny_manual <- synteny_manual %>%
  mutate(
    Start_1 = pmax(Start_1 - 400000, 1),
    End_1   = pmax(End_1 - 400000, 1),
    Start_2 = pmax(Start_2 - 400000, 1),
    End_2   = pmax(End_2 - 400000, 1)
  )

# 3. Mark block 2 as an inversion (swap Start_2/End_2 for row 2 only)
synteny_manual[2, c("Start_2", "End_2")] <- synteny_manual[2, c("End_2", "Start_2")]

print("Updated karyotype_manual:")
print(karyotype_manual)
print("Updated synteny_manual (2nd block inverted):")
print(synteny_manual)

##4. Generate ideogram SVG 
svgfile <- "manual_data_ideogram.svg"
ideogram(karyotype = karyotype_manual, synteny = synteny_manual, output = svgfile)

# 5. Convert to PNG to load in adobe illustrator
if (file.exists(svgfile)) {
  convertSVG(svgfile, file = "manual_data_ideogram", device = "png")
  svg_content <- rsvg::rsvg(svgfile)
  grid::grid.raster(svg_content)
} else {
  stop("Error: unable to generate SVG file, check data inputs.")
}
