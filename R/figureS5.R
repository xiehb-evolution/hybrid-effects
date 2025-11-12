# Custom theme function (enhanced for better aesthetics)
sci_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Text elements
      plot.title = element_text(size = base_size*1.3, hjust = 0.5, face = "bold", margin = margin(b = 15)),
      plot.subtitle = element_text(size = base_size*1.1, hjust = 0.5),
      plot.caption = element_text(size = base_size*0.8, hjust = 0),
      
      # Axis elements
      axis.title = element_text(size = base_size*1.1, face = "bold"),
      axis.text = element_text(size = base_size*0.9, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      
      # Legend elements
      legend.title = element_text(size = base_size*1.0, face = "bold"),
      legend.text = element_text(size = base_size*0.9),
      legend.key.size = unit(1.0, "lines"),
      legend.background = element_rect(fill = "white", color = "gray90"),
      legend.margin = margin(3, 3, 3, 3),
      
      # Panel elements
      panel.grid.major = element_line(color = "gray95", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),
      panel.background = element_rect(fill = "white"),
      
      # Facet elements
      strip.background = element_rect(fill = "#f5f5f5", color = "black", size = 0.5),
      strip.text = element_text(size = base_size*1.0, face = "bold", margin = margin(5, 0, 5, 0)),
      
      # Plot margins
      plot.margin = margin(8, 8, 8, 8)
    )
}

library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)  # Added cowplot for plot_grid function

# Data retrieval and processing
query <- paste0("SELECT chr, window, tag, lambda, WEIGHTED_FST, snps, ",
                "recombination_rate_all, mean_rec_all ",
                "FROM renew_complete_data_with_overall2 ",
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
# Define improved color mapping for tags (more visually distinct)
color_mapping <- c(
  "inbreeding depression" = "#3273dc", 
  "hybrid vigor" = "#23d160",          
  "hybrid depression" = "#ff3860"      
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
      x = "",  # Remove x-axis label as it will be shown only on the bottom plot
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
  # Stack the three plots vertically with 2:1:1 ratio and reduced spacing
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

# Remove individual legends from all plots
plot_chr1 <- plot_chr1 & theme(legend.position = "none")
plot_chr8 <- plot_chr8 & theme(legend.position = "none")
plot_chr6 <- plot_chr6 & theme(legend.position = "none")
plot_chr15 <- plot_chr15 & theme(legend.position = "none")

# Arrange in 2x2 grid with EFGH layout
# Using cowplot's plot_grid for precise control
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
# 添加共享图例
final_plot <- plot_grid(
  legend,
  combined_plot, 
  ncol = 1,
  rel_heights = c(0.05, 1)
)
final_plot
# Save the plot with appropriate dimensions for a 2x2 layout
# Using dimensions appropriate for a journal figure
ggsave("chromosome_EFGH_layout.png", final_plot, width = 10, height = 10, dpi = 300)
ggsave("chromosome_ABCD.pdf", final_plot, width = 12, height = 10)  
getwd()
# Display the final plot
final_plot












# Convert patchwork plots to grobs for use with cowplot
plot_chr1_grob <- patchwork::patchworkGrob(plot_chr1)
plot_chr3_grob <- patchwork::patchworkGrob(plot_chr3)

# Create the complete figure
# Assuming panel_1a through panel_3f are already defined elsewhere
figure1 <- plot_grid(
  panel_1a, panel_1b,
  panel_2a, panel_2b,
  panel_3d, panel_3e, panel_3f,
  plot_chr1_grob, plot_chr3_grob,  # Add chromosome plots as panels H and I
  ncol = 2,
  align = 'v',
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),  # Updated labels to include H and I
  label_size = 16,  # Adjust size as needed
  label_fontfamily = "sans", # Use a standard font
  label_fontface = "bold",
  hjust = -0.2,  # Adjust horizontal position of labels
  vjust = 1.1    # Adjust vertical position of labels
)
figure1
# Save the figure
ggsave("figure1_complete.pdf", figure1, width = 12, height = 15, units = "in")