# Shorter version
library(heatmaply)
library(plotly)
library(circlize)
heatmap_genes_plot <- function(dataframe, genes,titles) {
  # Extract the relevant columns based on the specified genes
  genes_data <- dataframe[, genes]
  
  # Calculate the correlation matrix
  correlation_matrix <- cor(genes_data, method = "pearson")
  
  library(ComplexHeatmap)
  col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  heatmap_result<-Heatmap(correlation_matrix,
        name = "Correlation",
        col = col_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 30),       # Adjust row names font size
        column_names_gp = gpar(fontsize = 30),    # Adjust column names font size
        column_title = titles,     # Add a title
        column_title_gp = gpar(fontsize = 35, fontface = "bold"),
          heatmap_legend_param = list(
          title_gp = gpar(fontsize = 0, fontface = "bold"),       # Legend title font size and style
          labels_gp = gpar(fontsize = 25)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = fill))
          grid.text(sprintf("%.2f", correlation_matrix[i, j]), x = x, y = y,
                    gp = gpar(fontsize = 20))
        })
  return(heatmap_result)
}



heatmap_genes_cor <- function(dataframe, genes) {
  # Extract the relevant columns based on the specified genes
  genes_data <- dataframe[, genes]
  
  # Calculate the correlation matrix
  correlation_matrix <- cor(genes_data, method = "pearson")
  return(correlation_matrix)
}



# Normalized
HPV_log2_cpm <- read.csv("all_batches_HPV_log2cpm_updated_HPVpos.csv")
norm_df <- HPV_log2_cpm

basic_genes <- c("e6","e7","e1", "e2", "e5", "l1", "l2")
#check all correlation
heatmap_genes_cor(norm_df, basic_genes)

norm_df$X
#visulize the heatmap
all_cohort<-heatmap_genes_plot(norm_df, basic_genes,"All cohorts HPV Genes")
all_cohort


png('Figure2E.png', units='in', width=10, height=10, res=300)
draw(all_cohort)
dev.off()
