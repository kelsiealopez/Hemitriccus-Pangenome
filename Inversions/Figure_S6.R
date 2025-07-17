if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

install.packages("DBI")
library(DBI)
BiocManager::install("clusterProfiler", force = TRUE)
# Load the package
library(clusterProfiler)



BiocManager::install("org.Hs.eg.db", force=TRUE)
library(org.Hs.eg.db)

gene_list <- c("ASNSD1",
               "ASDURF",
               "WDR75",
               "RAMP1",
               "RBM44",
               "LRRFIP1",
               "PRLH",
               "MLPH",
              "chiLan.LOC116789397",
              "SLC40A1",              
              "chiLan.LOC116789402",
              "chiLan.LOC116789705"  ,            
              "chiLan.LOC116789401",
              "UBE2F" ,             
              "RAB17",
              "MREG"
)


converted_genes <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


go_enrichment <- enrichGO(gene = converted_genes$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",  # Other options include "MF" and "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)


# View top results
head(summary(go_enrichment))

# Bar plot
barplot(go_enrichment, showCategory = 10)

# Dot plot
dotplot(go_enrichment, showCategory = 10)

dotplot(go_enrichment, showCategory = 20)




# Convert the enrichment result to a data frame
go_enrichment_df <- as.data.frame(go_enrichment)



library(ggplot2)

# Create a ggplot for the GO enrichment results
ggplot(go_enrichment_df, aes(x = reorder(Description, -Count), y = Count)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +  # Customize gradient colors
  coord_flip() +
  theme_classic() +
  labs(title = "GO Enrichment Analysis",
       y = "Gene Count",
       x = "GO Term",
       color = "Adjusted\np-value") +
  theme(axis.text.y = element_text(size = 10, hjust = 1))
