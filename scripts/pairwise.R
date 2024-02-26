# Load necessary libraries
library(data.table)
library(ggplot2)
library(tidyr)

# Function to compute pairwise L1 distances
compute_pairwise_L1 <- function(data) {
  distances <- as.matrix(dist(data, method = "manhattan"))
  return(distances)
}

# Read the matrix and take only the first 1000 rows
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
df <- fread(file_name, header = TRUE, sep = "\t")[1:2000, ]

# Select every second column starting from the fourth
df_selected <- df[, .SD, .SDcols = seq(4, ncol(df), by = 2)]

# Normalize rows function
normalize_rows <- function(df) {
  return(as.data.frame(t(apply(df, 1, scale))))
}

# Normalize the columns
normalize_columns <- function(df) {
  return(scale(df))
}
df_normalized <- as.data.frame(normalize_rows(df_selected))
#df_normalized <- df

# Calculate and plot pairwise L1 distances for each sequential sample
for (i in seq(1, ncol(df_normalized), by = 1)) {
  distances <- compute_pairwise_L1(df_normalized[, 1:i])
  distances_df <- as.data.frame(distances)
  distances_long <- pivot_longer(distances_df, cols = everything(), names_to = "pair", values_to = "distance")
  print(distances_long)

  # Plotting
  ggplot(distances_long, aes(x = distance)) +
    geom_histogram(bins = 100, fill = "blue", color = "black") +
    ggtitle(paste("Distribution of Pairwise L1 Distances for", i, "Columns")) +
    xlab("L1 Distance") +
    ylab("Frequency") +
    xlim(0,10)
  ggsave(paste("L1_Distances_", i, "_Columns.png", sep = ""))
}
