# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
# Read the data
file1 <- read_tsv(args[1])
file2 <- read_tsv(args[2])
head(file1)

# Select only the numeric columns
file1_numeric <- select(file1, where(is.numeric))
file2_numeric <- select(file2, where(is.numeric))

# Ensure both data frames have the same number of rows
if (nrow(file1_numeric) != nrow(file2_numeric)) {
    stop("The number of rows in the files does not match.")
}


# Select even columns. Adjust this for odd columns if needed
file1_even <- file1_numeric[2, seq(2, ncol(file1_numeric), by = 2)]
file2_even <- file2_numeric[2, seq(2, ncol(file2_numeric), by = 2)]

# Calculate Pearson correlations for even columns
correlations <- mapply(cor, as.data.frame(t(file1_even)), as.data.frame(t(file2_even)), MoreArgs = list(method = "pearson"))

# Create a data frame for ggplot
cor_df <- data.frame(correlation = correlations)
head(cor_df)
print(median(cor_df$correlation))

x11()

ggplot(cor_df, aes(x = factor(1), y = correlation)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, color = "blue", size = 2) +
    labs(title = "Pearson Correlations Between Corresponding Rows",
         x = "",
         y = "Correlation Coefficient") +
    theme_minimal()

Sys.sleep(Inf)
