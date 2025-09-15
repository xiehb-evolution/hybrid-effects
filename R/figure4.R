library(RMySQL)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales) 

# Define a clean, minimalist theme with white background and no grid
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

# Read data
data <- dbGetQuery(con, "SELECT * FROM stat_dev_from_sex_mean_in_fst_bin2")
data2 <- dbGetQuery(con, "SELECT * FROM hybrid_effect_analysis")
data3 <- dbGetQuery(con, "SELECT * FROM stat_dev_from_sex_mean_in_fst_bin")


# Plot A: Male and female differences across FST values
plotA <- ggplot(data) +
  geom_smooth(aes(x = fst_midpoint, y = malediff, color = "Male"), 
              method = "loess", span = 1.25, se = TRUE, linetype = "solid", size = 0.8) +
  geom_smooth(aes(x = fst_midpoint, y = femalediff, color = "Female"), 
              method = "loess", span = 1.25, se = TRUE, linetype = "solid", size = 0.8) +
  geom_point(aes(x = fst_midpoint, y = malediff, color = "Male"), 
             alpha = 0.8, size = 2, shape = 16) +
  geom_point(aes(x = fst_midpoint, y = femalediff, color = "Female"), 
             alpha = 0.8, size = 2, shape = 17) +
  scale_color_manual(values = c("Male" = "#3366CC", "Female" = "#CC3366"),
                     name = NULL) +
  labs(
    x = expression(F[ST]),
    y = "MPH (%)"
  ) +
  clean_theme +
  scale_y_continuous(labels = function(y) paste0(y * 100)) +
  theme(
    legend.position = c(0.2, 0.9),
    legend.background = element_rect(fill = "white", color = NA)
  )
plotA


# Plot B: Relationship between Lambda and sex differences
cor_result_B_male <- cor.test(data$malediff, data$lambda)
r_value_B_male <- cor_result_B_male$estimate
p_value_B_male <- cor_result_B_male$p.value
cor_result_B_female <- cor.test(data$femalediff, data$lambda)
r_value_B_female <- cor_result_B_female$estimate
p_value_B_female <- cor_result_B_female$p.value
data_long <- rbind(
  data.frame(diff = data$malediff, lambda = data$lambda, gender = "Male"),
  data.frame(diff = data$femalediff, lambda = data$lambda, gender = "Female")
)
print(range(data_long$lambda)) 
print(range(data_long$diff * 100)) 
plotB <- ggplot(data_long, aes(y = diff * 100, x = lambda, color = gender)) +
  geom_smooth(method = "lm", linetype = "solid", size = 0.8, se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Female" = "#CC3366", "Male" = "#3366CC")) +
  geom_point(alpha = 0.8, size = 1.5) +
  labs(y = "MPH (%)", x = expression(lambda), color = "Sex") +
  annotate(
    "text", 
    y = 2.61,
    x = 33.3,
    label = sprintf("Male: r = %.3f, p = %.2e\nFemale: r = %.3f, p = %.2e", 
                    r_value_B_male, p_value_B_male, r_value_B_female, p_value_B_female),
    size = 2.5, 
    hjust = 0
  ) +
  coord_cartesian( 
    xlim = range(data_long$lambda, na.rm = TRUE),
    ylim = range(data_long$diff * 100, na.rm = TRUE)
  ) +
  clean_theme +
  theme(legend.position = c(0.8, 0.9))
plotB


# Plot C: Heterozygote advantage ratio vs sex differences
merged_data <- data %>%
  inner_join(data2, by = "fst_bin")
cor_result_C_male <- cor.test(merged_data$het_advantage_ratio, merged_data$malediff)
r_value_C_male <- cor_result_C_male$estimate
p_value_C_male <- cor_result_C_male$p.value
cor_result_C_female <- cor.test(merged_data$het_advantage_ratio, merged_data$femalediff)
r_value_C_female <- cor_result_C_female$estimate
p_value_C_female <- cor_result_C_female$p.value

# Create a long-format dataset for plotting
data_long <- rbind(
  data.frame(het_advantage_ratio = merged_data$het_advantage_ratio, 
             diff = merged_data$malediff, 
             gender = "Male"),
  data.frame(het_advantage_ratio = merged_data$het_advantage_ratio, 
             diff = merged_data$femalediff, 
             gender = "Female")
)

plotC <- ggplot(data_long, aes(x = het_advantage_ratio, y = diff*100, color = gender)) +
  geom_smooth(method = "lm", linetype = "solid", size = 0.8, se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Female" = "#CC3366", "Male" = "#3366CC")) +
  geom_point(alpha = 0.8, size = 1.5) +
  labs(
    x = "Heterozygote advantage/disadvantage ratio", 
    y = "MPH (%)", 
    color = "Sex"
  ) +
  annotate(
    "text", 
    x = 1,
    y = 2.6,
    label = sprintf("Male: r = %.3f, p = %.2e\nFemale: r = %.3f, p = %.2e", 
                    r_value_C_male, p_value_C_male, r_value_C_female, p_value_C_female),
    size = 2.5, hjust = 0
  ) +
  clean_theme +
  theme(legend.position = c(0.8, 0.9))
plotC




# Plot D: Male and female phenotypic differences
sql = "select fst_midpoint, malediff as diff, 'Male' as sex from stat_dev_from_sex_mean_in_fst_bin 
       union select fst_midpoint, femalediff as diff, 'Female' as sex from stat_dev_from_sex_mean_in_fst_bin"
figure4ddata = dbGetQuery(con, sql)
plotD <- ggplot(data = figure4ddata, aes(x = fst_midpoint, y = diff*100, colour = sex)) + 
  scale_color_manual(values = c("Female" = "#CC3366", "Male" = "#3366CC")) +
  geom_smooth(method = "loess", span = 1.25) +
  geom_point(size = 1.5) + 
  labs(
    y = "Homozygote phenotypic difference (%)",
    x = expression(F[ST]),
    colour = "Sex" 
  ) +
  clean_theme + 
  theme(legend.position = c(0.2, 0.85))
plotD


add_outside_labels <- function(plot, label, 
                               x_offset = 0.05,  
                               y_offset = 0.01,  
                               plot_width = 0.95, 
                               label_size = 11,
                               left_margin = 10, 
                               top_margin = 10,  
                               bottom_margin = 5, 
                               right_margin = 5) {

  plot_adjusted <- plot + 
    theme(
      plot.margin = margin(top_margin, right_margin, bottom_margin, left_margin),
      axis.title.y = element_text(margin = margin(r = 10))
    )

  plot_with_label <- ggdraw() + 
    draw_plot(plot_adjusted, x = x_offset, y = y_offset, 
              width = plot_width, height = 0.98 - y_offset) +
    draw_label(label, x = x_offset/3, y = 0.99, 
               hjust = 0, vjust = 1, 
               size = label_size, fontface = "bold")
  
  return(plot_with_label)
}
plotA_labeled <- add_outside_labels(plotA, "A")
plotB_labeled <- add_outside_labels(plotB, "B")
plotC_labeled <- add_outside_labels(plotC, "C")
plotD_labeled <- add_outside_labels(plotD, "D")
combined_plot <- plot_grid(
  plotA_labeled, plotB_labeled, 
  plotC_labeled, plotD_labeled,
  ncol = 2,
  align = "hv"  
)

plotA
plotB
plotC
plotD
combined_plot
ggsave("combined_figure.png", combined_plot, width = 10, height = 6, dpi = 300)
ggsave("figure4.pdf", combined_plot, width = 10, height = 6, device = cairo_pdf)
