# Load necessary libraries
library(data.table)
library(Rtsne)
library(ggplot2)

# Normalize rows function
normalize_rows <- function(df) {
  return(as.data.frame(t(apply(df, 1, scale))))
}

# Normalize the columns
normalize_columns <- function(df) {
  return(scale(df))
}


# Read the matrix and take only the first 1000 rows
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
df <- fread(file_name, header = TRUE, sep = "\t")[1:3000, ]

# Select every second column starting from the fourth
df_selected <- df[, .SD, .SDcols = seq(4, ncol(df), by = 2)]

# Normalize the rows
df_normalized <- normalize_columns(df_selected)
#df_normalized <- df

# Perform and plot t-SNE for each sequential sample
for (i in seq(1, ncol(df_normalized), by = 1)) {
  # Remove duplicate rows
  df_unique <- unique(df_normalized[, 1:i])
  print(df_unique[0:10])

  # Run t-SNE only if there are enough unique rows
  #if (nrow(df_unique) > 1) {
    set.seed(42)  # For reproducibility
    tsne_result <- Rtsne(df_unique, dims = 2, perplexity = 30, verbose = FALSE)

    # Prepare data for ggplot
    tsne_data <- as.data.frame(tsne_result$Y)
    colnames(tsne_data) <- c("TSNE1", "TSNE2")
    
    # Plotting
    ggplot(tsne_data, aes(x = TSNE1, y = TSNE2)) +
      geom_point() +
      ggtitle(paste("t-SNE Visualization for", i, "Columns")) +
      xlab("t-SNE Dimension 1") +
      ylab("t-SNE Dimension 2")
    ggsave(paste("TSNE_", i, "_Columns.png", sep = ""))
  }
}
