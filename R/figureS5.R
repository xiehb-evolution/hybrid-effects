# Custom theme function (enhanced for better aesthetics)
sci_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size*1.3, hjust = 0.5, face = "bold", margin = margin(b = 15)),
      plot.subtitle = element_text(size = base_size*1.1, hjust = 0.5),
      plot.caption = element_text(size = base_size*0.8, hjust = 0),
      
      axis.title = element_text(size = base_size*1.1, face = "bold"),
      axis.text = element_text(size = base_size*0.9, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      
      legend.title = element_text(size = base_size*1.0, face = "bold"),
      legend.text = element_text(size = base_size*0.9),
      legend.key.size = unit(1.0, "lines"),
      legend.background = element_rect(fill = "white", color = "gray90"),
      legend.margin = margin(3, 3, 3, 3),
      
      panel.grid.major = element_line(color = "gray95", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),
      panel.background = element_rect(fill = "white"),
      
      strip.background = element_rect(fill = "#f5f5f5", color = "black", size = 0.5),
      strip.text = element_text(size = base_size*1.0, face = "bold", margin = margin(5, 0, 5, 0)),
      
      plot.margin = margin(8, 8, 8, 8)
    )
}

library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
# Data retrieval and processing
query <- paste0("SELECT chr, window, tag, lambda, WEIGHTED_FST, snps, ",
                "recombination_rate_all, mean_rec_all ",
                "FROM renew_complete_data_with_overall ",
                "WHERE Sex = 'Overall'")
hybrid_effect_data <- dbGetQuery(con, query)
hybrid_effect_data$chr <- as.numeric(hybrid_effect_data$chr)
hybrid_effect_data$tag <- factor(hybrid_effect_data$tag, 
                                 levels = c("inbreeding depression", "hybrid vigor", "hybrid depression"))
chr_data <- hybrid_effect_data %>% 
  filter(chr >= 1, chr <= 18) %>% 
  mutate(
    chr = factor(chr),
    position_mb = (window - 1) * 0.1,
    tag = factor(tag, levels = c("inbreeding depression", "hybrid vigor", "hybrid depression"))
  )
color_mapping <- c(
  "inbreeding depression" = "#E69F00",
  "hybrid vigor" = "#009E73",          
  "hybrid depression" = "#0072B2"      
)
# Define shape mapping
shape_mapping <- c(
  "inbreeding depression" = 21, 
  "hybrid vigor" = 22,          
  "hybrid depression" = 24      
)

# Create functions to generate three-panel plot for each chromosome
create_chr_plots <- function(chr_num) {
  # Filter data for the specific chromosome
  chr_specific_data <- filter(chr_data, chr == chr_num)
  # Plot 1: Enhanced scatter plot with points colored by tag
  p1 <- ggplot(chr_specific_data, 
               aes(x = position_mb, y = lambda)) +
    geom_point(aes(shape = tag, fill = tag, color = tag),
               size = 2, stroke = 0.5, alpha = 0.8) +
    scale_shape_manual(values = shape_mapping) +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    labs(
      x = "", 
      y = expression(lambda)
    ) +
    sci_theme() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank()  
    )
  
  # Plot 2: Enhanced line plot for WEIGHTED_FST
  p2 <- ggplot(chr_specific_data, 
               aes(x = position_mb, y = WEIGHTED_FST)) +
    geom_line(color = "#1a237e", size = 1.0) + 
    geom_area(fill = "#7986cb", alpha = 0.2) + 
    labs(
      x = "",  
      y = "FST"
    ) +
    sci_theme() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank() 
    )
  
  # Plot 3: Bar chart for recombination rate
  p3 <- ggplot(chr_specific_data, 
               aes(x = position_mb, y = recombination_rate_all)) +
    geom_line(color = "darkred", size = 0.8) +
    labs(
      x = "Chromosomal position (Mb)", 
      y = "Recombination rate"
    ) +
    sci_theme() +
    theme(
      legend.position = "none"
    )
  stacked_plot <- p1 / p2 / p3 +
    plot_layout(heights = c(2, 1, 1), guides = "collect") +
    plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) &
    theme(
      plot.margin = margin(1, 5, 1, 5),  
      legend.position = "top",           
      legend.justification = "center",   
      legend.box.margin = margin(0, 0, 5, 0),  
      legend.title = element_blank()  
    )
  
  return(stacked_plot)
}

# Create plots for chromosome 1 and 3
# Modify sci_theme temporarily to include legend position
chr_theme <- function(base_size = 12) {
  sci_theme(base_size) + 
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
}

plot_chr1 <- create_chr_plots(1) + 
  plot_annotation(title = "Chromosome 1", theme = chr_theme())
plot_chr1

plot_chr8 <- create_chr_plots(8) + 
  plot_annotation(title = "Chromosome 8", theme = chr_theme())
plot_chr8

plot_chr6 <- create_chr_plots(6) + 
  plot_annotation(title = "Chromosome 6", theme = chr_theme())
plot_chr6

plot_chr15 <- create_chr_plots(15) + 
  plot_annotation(title = "Chromosome 15", theme = chr_theme())
plot_chr15


# Building on your existing code to create an EFGH layout grid
# Create the base chromosome plots without titles
plot_chr1 <- create_chr_plots(1)
plot_chr8 <- create_chr_plots(8)
plot_chr6 <- create_chr_plots(6)
plot_chr15 <- create_chr_plots(15)

# Create a shared legend for all plots to avoid repetition
# Extract legend from one of the plots
legend <- get_legend(plot_chr1)
plot_chr1 <- plot_chr1 & theme(legend.position = "none")
plot_chr8 <- plot_chr8 & theme(legend.position = "none")
plot_chr6 <- plot_chr6 & theme(legend.position = "none")
plot_chr15 <- plot_chr15 & theme(legend.position = "none")

combined_plot <- plot_grid(
  plot_chr1, plot_chr8,
  plot_chr6, plot_chr15,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 16,
  label_fontface = "bold",
  align = "hv",
  axis = "tblr",
  hjust = 0.8,
  vjust = 1.5,
  margin = unit(0.5, "cm")
) + theme(
  plot.margin = margin(10, 10, 10, 10)
)
final_plot <- plot_grid(
  legend,
  combined_plot, 
  ncol = 1,
  rel_heights = c(0.05, 1)
)
final_plot





