##############################################################################

# Publication-Ready Analysis of Hybrid Effects in Sexual Reproduction

##############################################################################

# Required packages
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(viridis)
library(grid)
library(cowplot)
library(scales)
library(RColorBrewer)
library(DescTools)
library(fitdistrplus)
library(RMySQL)
library(sqldf)

# Set global options
options(scipen = 999) 
set.seed(42)
sci_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size*1.2, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size*1.1, hjust = 0.5),
      plot.caption = element_text(size = base_size*0.8, hjust = 0),
      
      axis.title = element_text(size = base_size*1.1, face = "bold"),
      axis.text = element_text(size = base_size*0.9, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),  
      axis.ticks = element_line(color = "black", size = 0.5),
      
      legend.title = element_text(size = base_size*1.0, face = "bold"),
      legend.text = element_text(size = base_size*0.9),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "white", color = "gray90"),
      legend.margin = margin(2, 2, 2, 2),
      
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),  
      
      strip.background = element_rect(fill = "gray95", color = "black", size = 0.5),
      strip.text = element_text(size = base_size*1.0, face = "bold"),
      
      plot.margin = margin(5, 5, 5, 5)
    )
}

# Function to format p-values scientifically
format_pvalue <- function(p, digits = 2) {
  if (p < 0.001) {
    return(paste0("p < 0.001"))
  } else {
    return(paste0("p = ", format(p, digits = digits)))
  }
}

# Function to add panel labels
add_panel_labels <- function(plot, label, x = -0.07, y = 1.07, size = 12, fontface = "bold") {
  plot + annotation_custom(
    grob = textGrob(label, x = unit(x, "npc"), y = unit(y, "npc"), 
                    gp = gpar(fontsize = size, fontface = fontface))
  )
}

# Figure SNPs count by hybrid effect type with statistical comparison
# Get data for hybrid effect types (inbreeding depression, hybrid vigor, hybrid depression)
query <- paste0("SELECT chr, window, tag, lambda, WEIGHTED_FST, snps, ",
                "recombination_rate_all, mean_rec_rate ",
                "FROM renew_complete_data_with_overall ",
                "WHERE Sex = 'Overall'")
hybrid_effect_data <- dbGetQuery(con, query)
# Data preprocessing
hybrid_effect_data$chr <- as.numeric(hybrid_effect_data$chr)
hybrid_effect_data$tag <- factor(hybrid_effect_data$tag, 
                                 levels = c("inbreeding depression", "hybrid vigor", "hybrid depression"))

panel_f_data <- hybrid_effect_data
panel_f_data$snps_winsorized <- Winsorize(panel_f_data$snps, probs = c(0, 0.95))
panel_3f <- ggplot(panel_f_data, aes(x = tag, y = snps_winsorized, fill = tag)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, outlier.alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "SNPs count (winsorized at 95%)",
    title = "SNP density across hybrid effect types"
  ) +
  scale_x_discrete(labels = c("Inbreeding\ndepression", "Hybrid\nvigor", "Hybrid\ndepression")) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1, color = "black") +
  stat_compare_means(comparisons = list(
    c("inbreeding depression", "hybrid vigor"),
    c("hybrid vigor", "hybrid depression"),
    c("inbreeding depression", "hybrid depression")),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    size = 3) +
  sci_theme() +
  theme(legend.position = "none")
panel_3f

