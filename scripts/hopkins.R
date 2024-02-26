# Load necessary libraries
library(pracma)
library(data.table)
library(hopkins)

# Read the matrix
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
df <- fread(file_name, header = TRUE, sep = "\t")[1:8000,]
#head(df)
#df <- df[, .SD, .SDcols = -1]
df <- df[, .SD, .SDcols = seq(4, ncol(df), by = 2)]

df_normalized <- df
head(df)

# Sequential sampling for Hopkins statistic
for (i in seq(1, ncol(df_normalized), by = 1)) {
    print(hopkins(df_normalized[, 1:i]))
}

# Plotting
#plot(seq(4, ncol(df), by = 2), hopkins_values, type = "b", xlab = "Column Index", ylab = "Hopkins Statistic", main = "Hopkins Statistic as a Function of Column Selection")

