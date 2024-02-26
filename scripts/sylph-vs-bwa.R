library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggsci)

args <- commandArgs(trailingOnly = TRUE)

plot_nature_style <- function(data, title = "Your Plot Title Here") {

  # Ensure the variables are interpreted as symbols for ggplot

    head(data)
    data = filter(data, contigLen.x < 5500)
    p <- ggplot(data = data, aes(y = QC.reads.DRR310876_R1.fastq.gz, x = scaffolds.fasta.DRR310876_R1.fastq.gz.bam)) +
    geom_point(size=0.5, alpha=0.6, color="black") +
   # geom_smooth(method = "lm", color="red") +
    scale_x_log10() +
    scale_y_log10() +
    stat_cor(method = "spearman", label.x = 0, label.y = 2) +
    theme_bw()
    #theme(panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.border = element_blank(),
    #panel.background = element_blank()) 

    
    #theme_hc()
    #theme_minimal(base_size = 14) +
    #theme(text = element_text(family = "Times"),
    #      axis.title = element_text(face = "bold"),
    #      panel.grid.major = element_blank(),
    #      panel.grid.minor = element_blank(),
    #      panel.background = element_rect(fill = "white", colour = "black")) +
    #stat_cor(method = "pearson", label.x = 0, label.y = 3) +
    #labs(x = xlab, y = ylab, title = title, subtitle = "Linear Regression and Correlation Coefficients")

  # Return the plot
  return(p)
}

df1 <- read.table(args[1], sep='\t', header=TRUE)
df2 <- read.table(args[2], sep='\t', header=TRUE)
#df2 <- rename(df2, Contig_name = contigName)
df <- left_join(df1,df2,by = 'contigName')

head(df)

x11()
p = plot_nature_style(df)
p + xlab("BWA coverage") + 
    ylab("Sylph-bincov coverage")
Sys.sleep(Inf)


