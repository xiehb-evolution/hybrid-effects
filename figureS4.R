library(RMySQL) 
library(dplyr) 
library(ggplot2) 
library(gridExtra)
library(cowplot)
library(scales)
library(extrafont)
#font_import()
fonts()
clean_theme <- theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 20)), 
    axis.text = element_text(size = 8, color = "black"),
    axis.text.y = element_text(margin = margin(r = 10)), 
    axis.line = element_line(color = "black", size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "lines"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 10, 5, 30) 
  )

query <- "SELECT chr, window, WEIGHTED_FST as fst FROM fst100k WHERE chr BETWEEN 1 AND 18"
fst_data <- dbGetQuery(con, query)
fst_data$fst[fst_data$fst < 0] <- 0
str(fst_data)
fst_summary <- summary(fst_data$fst)
print(fst_summary)
chromosome_means <- aggregate(fst ~ chr, fst_data, mean)
print(chromosome_means)
########################################################################################
median_fst <- median(fst_data$fst)
dens <- density(fst_data$fst)
max_density <- max(dens$y)
y_max <- max_density * 1.15
pdensity <- ggplot(fst_data, aes(x = fst)) +
  
  geom_density(fill = "#8ECAE6", 
               color = "#219EBC",
               alpha = 0.85,
               linewidth = 0.8) +
  
  geom_vline(
    xintercept = median_fst,
    linetype = "dashed", 
    color = "#FB8500",
    linewidth = 1.0
  ) +
  
  annotate(
    "text",
    x = median_fst * 1.15,
    y = max_density * 0.92,
    label = paste0("Median = ", round(median_fst, 4)),
    color = "#FB8500",
    size = 4.2,
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_x_continuous(
    limits = c(0, max(fst_data$fst)),
    expand = c(0, 0),
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::comma_format(accuracy = 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = c(0, 0)
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5), 
    axis.ticks = element_line(color = "black", linewidth = 0.5), 
    
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(color = "grey40"),
    plot.title = element_text(
      face = "bold", 
      hjust = 0.5, 
      size = 14,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "grey50",
      margin = margin(b = 12)
    ),
    plot.margin = margin(10, 15, 10, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    x = expression(F[ST]),
    y = "Density"
  )
print(pdensity)
pdensity <- pdensity + theme(plot.margin = margin(t = 15, r = 0, b = 0, l = 10, unit = "pt"))
figureS2 <- plot_grid(
  pdensity,
  ncol = 1,  
  align = 'hv',
  label_size = 16,
  label_fontfamily = "sans",
  label_fontface = "bold",
  hjust = 0,
  vjust = 1.2
)
figureS2
