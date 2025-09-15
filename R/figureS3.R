##############################################################################
# Publication-Ready Analysis of Hybrid Effects in Sexual Reproduction
# Figures for Scientific Journal Article
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
options(scipen = 999)
set.seed(42) 
sci_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size*1.2, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size*1.1, hjust = 0.5),

      axis.title = element_text(size = base_size*1.1, face = "bold"),
      axis.text = element_text(size = base_size*0.9, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),  # 控制坐标轴线
      axis.ticks = element_line(color = "black", size = 0.5),

      legend.title = element_text(size = base_size*1.0, face = "bold"),
      legend.text = element_text(size = base_size*0.9),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "white", color = "gray90"),
      legend.margin = margin(2, 2, 2, 2),
      
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),  # 移除整个面板边框
      
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

##############################################################################
# FIGURE S1: Hybrid Effect Size Distributions
##############################################################################
# Panel A: Male hybrid effect size distribution
sql = "SELECT abs(malemutantdev-maledev) as effect FROM window100k_single_site_trait_stat_mutant_deviation_from_mean where chr<23 having effect<=0.3"
male_effect_data = dbGetQuery(con, sql)
male_fit = fitdistr(male_effect_data$effect, "exponential") 
male_lambda = male_fit$estimate

panel_1a = ggplot(data=male_effect_data, aes(x=effect)) + 
  geom_density(color="#3366CC", fill="#3366CC", alpha=0.2, linewidth=1) + 
  xlim(0, 0.3) + 
  sci_theme() +
  labs(
    x = "Hybrid effect size", 
    y = "Density",
    title = "Male hybrid effect distribution"
  ) +
  annotate("text", x = 0.22, y = max(density(male_effect_data$effect)$y)*0.9, 
           label = paste0("λ = ", round(male_lambda, 2)), 
           hjust = 0, size = 3.5)
panel_1a

# Panel B: Female hybrid effect size distribution
sql = "SELECT abs(femalemutantdev-femaledev) as effect FROM `window_single_site_trait_stat_mutant_deviation_from_mean` where chr<23 having effect<=0.3"
female_effect_data = dbGetQuery(con, sql)
female_fit = fitdistr(female_effect_data$effect, "exponential") 
female_lambda = female_fit$estimate

panel_1b = ggplot(data=female_effect_data, aes(x=effect)) + 
  geom_density(color="#CC6677", fill="#CC6677", alpha=0.2, linewidth=1) + 
  xlim(0, 0.3) + 
  sci_theme() +
  labs(
    x = "Hybrid effect size", 
    y = "Density",
    title = "Female hybrid effect distribution"
  ) +
  annotate("text", x = 0.22, y = max(density(female_effect_data$effect)$y)*0.9, 
           label = paste0("λ = ", round(female_lambda, 2)), 
           hjust = 0, size = 3.5)

panel_1b

panel_1a <- panel_1a + theme(plot.title = element_blank())
panel_1b <- panel_1b + theme(plot.title = element_blank())
panel_1a <- panel_1a + theme(plot.margin = margin(t = 15, r = 0, b = 0, l = 10, unit = "pt"))
panel_1b <- panel_1b + theme(plot.margin = margin(t = 15, r = 0, b = 0, l = 10, unit = "pt"))
figureS1 <- plot_grid(
  panel_1a, panel_1b,
  ncol = 2,  
  align = 'hv',
  labels = c("A", "B"),
  label_size = 16,
  label_fontfamily = "sans",
  label_fontface = "bold",
  hjust = 0,
  vjust = 1.2
)
figureS1
#ggsave("figureS1.pdf", figureS1, width = 12, height = 5, device = cairo_pdf)
getwd()

