library(ggplot2)
library(RMySQL)
library(sqldf)
library(dplyr)
library(gridExtra)
library(tidyr)
library(scales)
library(viridis)
library(cowplot)
sci_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size*1.2, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size*1.1, hjust = 0.5),
      
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
zero_grey <- "grey80"
add_panel_labels <- function(plot, label, x = -0.07, y = 1.07, size = 12, fontface = "bold") {
  plot + annotation_custom(
    grob = textGrob(label, x = unit(x, "npc"), y = unit(y, "npc"), 
                    gp = gpar(fontsize = size, fontface = fontface))
  )
}

format_pvalue <- function(p, digits = 2) {
  if (p < 0.001) {
    return(paste0("p < 0.001"))
  } else {
    return(paste0("p = ", format(p, digits = digits)))
  }
}


#######################################################################################################
#Separate statistics for males and females
lambda_data <- dbGetQuery(con, "
    SELECT a.chr, a.window, a.fst_bin, a.fst_start, a.fst_end, a.WEIGHTED_FST, 
           a.lambda as lambda_male, 
           b.lambda as lambda_female
    FROM renew_complete_data_with_overall a 
    JOIN renew_complete_data_with_overall b 
      ON a.chr = b.chr AND a.window = b.window 
    WHERE a.sex = 'Male' AND b.sex = 'Female' AND a.chr < 23
") %>%
  mutate(fst_midpoint = (fst_start + fst_end) / 2) %>%
  group_by(fst_midpoint) %>%
  slice(1) %>%  
  ungroup()

lambda_summary_bin <- lambda_data %>%
  group_by(fst_bin) %>%
  summarise(
    across(c(lambda_male, lambda_female), 
           list(mean = mean, se = ~sd(.x)/sqrt(n())), 
           .names = "{.fn}_{.col}"),
    mean_fst_midpoint = mean(fst_midpoint),
    count = n()
  ) %>%
  arrange(mean_fst_midpoint)

plot_A <- ggplot(lambda_data, aes(x = fst_midpoint)) +
  geom_smooth(aes(y = lambda_male, color = "Male", fill = "Male"),
              method = "loess", span = 0.9, se = TRUE, alpha = 0.1, size = 1.0) +
  geom_smooth(aes(y = lambda_female, color = "Female", fill = "Female"),
              method = "loess", span = 0.9, se = TRUE, alpha = 0.1, size = 1.0) +
 
  geom_point(aes(y = lambda_male, color = "Male", shape = "Male"), 
             alpha = 0.7, size = 3.5) +
  geom_point(aes(y = lambda_female, color = "Female", shape = "Female"), 
             alpha = 0.7, size = 3.5) +
  
  scale_color_manual(name = "Sex",
                     values = c(Male = "#4682B4", Female = "#FF69B4"),
                     guide = guide_legend(override.aes = list(
                       shape = c(16, 17),
                       linetype = c(1, 1)
                     ))) +
  scale_fill_manual(name = "Sex", 
                    values = c(Male = "#4682B4", Female = "#FF69B4"),
                    guide = "none") +
  scale_shape_manual(values = c(Male = 16, Female = 17), guide = "none") +
  scale_linetype_manual(values = c(Male = 1, Female = 1), guide = "none") +
  
  labs(x = expression(F[ST]), y = expression(lambda)) +
  sci_theme(11) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.5),
    legend.key = element_blank()
  )
plot_A



#######################################################################################################
# Get rare SNP data (FST grouping)
sql_rare <- "SELECT round(a.WEIGHTED_FST/0.01,0) as fst,
2 * sum(((b.malehet1)/(b.malehet1+b.malehomo1) - (b.femalehet1)/(b.femalehet1+b.femalehomo1))/((b.malehet1)/(b.malehet1+b.malehomo1) + (b.femalehet1)/(b.femalehet1+b.femalehomo1)))/count(b.bias) as bias,
sum(case b.bias>0 when 1 then 1 end)/count(*) as male_bias_ratio,
sum(a.snps) as snps1
FROM renew_complete_data_with_overall a, Rare_mutation_sex_bias b
where a.chr=b.chr and a.window=b.window and a.sex='Overall' and a.weighted_fst<=0.5
group by round(a.WEIGHTED_FST/0.01,0)"
rare_result <- dbGetQuery(con, sql_rare)

# Get common SNP data (FST grouping)
sql_common <- "SELECT round(a.WEIGHTED_FST/0.01,0) as fst,
2 * sum(((b.malehet2)/(b.malehet2+b.malehomo2) - (b.femalehet2)/(b.femalehet2+b.femalehomo2))/((b.malehet2)/(b.malehet2+b.malehomo2) + (b.femalehet2)/(b.femalehet2+b.femalehomo2)))/count(b.bias) as bias,
sum(case b.bias>0 when 1 then 1 end)/count(*) as male_bias_ratio,
sum(a.snps) as snps1
FROM renew_complete_data_with_overall a, Rare_mutation_sex_bias b
where a.chr=b.chr and a.window=b.window and a.sex='Overall' and a.weighted_fst<=0.5
group by round(a.WEIGHTED_FST/0.01,0)"
common_result <- dbGetQuery(con, sql_common)

# Figure D: Relationship between female heterozygous loss of rare SNPs and FST
rare_result$fst_value <- rare_result$fst * 0.01 
plot_D <- ggplot(data=rare_result, aes(x=fst_value, y=bias)) +
  geom_smooth(method="lm", color="#e74c3c", se=TRUE, size=0.8) +
  labs(
    x = expression(F[ST]), 
    y = "Female heterozygote deficiency"
  ) +
  geom_point(size=3.5, alpha=0.8, color="#3498db") +
  sci_theme(11)
plot_D

# Calculating correlations for rare SNPs
rare_cor <- cor.test(rare_result$fst, rare_result$bias)
plot_D <- plot_D + 
  annotate("text", x = min(rare_result$fst_value), y = max(rare_result$bias),
           label = sprintf("r = %.3f, %s", rare_cor$estimate, format_pvalue(rare_cor$p.value)),
           hjust = 0, vjust = 1, size = 3.5)
plot_D


# Figure E: Relationship between gender differences in common SNPs and FST
common_result$fst_value <- common_result$fst * 0.01  
plot_E <- ggplot(data=common_result, aes(x=fst_value, y=bias)) +
  geom_smooth(method="lm", color="#9b59b6", se=TRUE, size=0.8) +
  labs(
    x = expression(F[ST]), 
    y = "Female heterozygote deficiency"
  ) +
  geom_point(size=3.5, alpha=0.8, color="#2ecc71") +
  sci_theme(11)

# Calculating correlations for common SNPs
common_cor <- cor.test(common_result$fst, common_result$bias)
plot_E <- plot_E + 
  annotate("text", x = min(common_result$fst_value), y = max(common_result$bias),
           label = sprintf("r = %.3f, %s", common_cor$estimate, format_pvalue(common_cor$p.value)),
           hjust = 0, vjust = 1, size = 3.5)
plot_E

# Get raw data of rare SNPs (chromosome-wide distribution)
sql_rare_raw <- "SELECT a.chr, a.window, a.WEIGHTED_FST as fst,
                ((b.malehet1)/(b.malehet1+b.malehomo1) - (b.femalehet1)/(b.femalehet1+b.femalehomo1)) as raw_diff,
                2 * ((b.malehet1)/(b.malehet1+b.malehomo1) - (b.femalehet1)/(b.femalehet1+b.femalehomo1)) / 
                ((b.malehet1)/(b.malehet1+b.malehomo1) + (b.femalehet1)/(b.femalehet1+b.femalehomo1)) as bias
                FROM renew_complete_data_with_overall a, Rare_mutation_sex_bias b
                WHERE a.chr=b.chr AND a.window=b.window AND a.sex='Overall' AND a.chr BETWEEN 1 AND 18"
rare_raw <- dbGetQuery(con, sql_rare_raw)

# Get raw data of common SNPs (chromosome-wide distribution)
sql_common_raw <- "SELECT a.chr, a.window, a.WEIGHTED_FST as fst,
                  ((b.malehet2)/(b.malehet2+b.malehomo2) - (b.femalehet2)/(b.femalehet2+b.femalehomo2)) as raw_diff,
                  2 * ((b.malehet2)/(b.malehet2+b.malehomo2) - (b.femalehet2)/(b.femalehet2+b.femalehomo2)) / 
                  ((b.malehet2)/(b.malehet2+b.malehomo2) + (b.femalehet2)/(b.femalehet2+b.femalehomo2)) as bias
                  FROM renew_complete_data_with_overall a, Rare_mutation_sex_bias b
                  WHERE a.chr=b.chr AND a.window=b.window AND a.sex='Overall' AND a.chr BETWEEN 1 AND 18"
common_raw <- dbGetQuery(con, sql_common_raw)

rare_raw$chr <- factor(rare_raw$chr, levels = 1:18)
common_raw$chr <- factor(common_raw$chr, levels = 1:18)
rare_raw$bias[is.na(rare_raw$bias)] <- 0
common_raw$bias[is.na(common_raw$bias)] <- 0

rare_raw <- rare_raw %>%
  group_by(chr) %>%
  arrange(window) %>%
  mutate(pos_index = row_number()) %>%
  ungroup()

common_raw <- common_raw %>%
  group_by(chr) %>%
  arrange(window) %>%
  mutate(pos_index = row_number()) %>%
  ungroup()


#######################################################################################################
# Calculate bias range for rare SNP and common SNP data
all_bias_values <- c(rare_raw$bias, common_raw$bias)
bias_range <- range(all_bias_values, na.rm = TRUE)
bias_limits <- c(-max(abs(bias_range)), max(abs(bias_range)))

female_color <- "#FF5252" 
male_color <- "#4169E1" 
zero_grey <- "#F5F5F5" 

# Figure B: Rare SNP Heatmap - Uses the same legend criteria
plot_B <- ggplot(rare_raw, aes(x = pos_index/10, y = chr, fill = bias)) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    low = female_color,
    mid = zero_grey,
    high = male_color,
    limits = bias_limits,
    name = "Sex bias",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.5,
      ticks.linewidth = 0.5,
      barwidth = 0.6,
      barheight = 6
    )
  ) +
  labs(
    x = "Position (Mb)",
    y = "Chromosome",
    title = "Rare SNPs"
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  sci_theme(12) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    legend.position = c(0.97, 0.4),
    legend.justification = c(1, 0.5),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.box.just = "right",
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),
    legend.title.align = 0.5,
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
plot_B

# Figure C: Common SNP heat map
plot_C <- ggplot(common_raw, aes(x = pos_index/10, y = chr, fill = bias)) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    low = female_color,
    mid = zero_grey,
    high = male_color,
    limits = bias_limits,
    name = "Sex bias",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.5,
      ticks.linewidth = 0.5,
      barwidth = 0.6,
      barheight = 6
    )
  ) +
  labs(
    x = "Position (Mb)",
    y = "Chromosome",
    title = "Common SNPs"
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  sci_theme(12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    legend.position = c(0.97, 0.4),
    legend.justification = c(1, 0.5),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.box.just = "right",
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),
    legend.title.align = 0.5,
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
plot_C



plot_A <- add_panel_labels(plot_A, "A")
plot_B <- add_panel_labels(plot_B, "B")
plot_C <- add_panel_labels(plot_C, "C")
plot_D <- add_panel_labels(plot_D, "D")
plot_E <- add_panel_labels(plot_E, "E")
print(plot_A)
print(plot_B)
print(plot_C)
print(plot_D)
print(plot_E)
library(cowplot)

figure2 <- plot_grid(
  plot_A, NULL, plot_B,
  plot_C, plot_D,
  plot_E,
  ncol = 3,
  align = 'hv',
  labels = c("A", "B", "C", "D", "E", "F"), 
  label_size = 16,
  label_fontfamily = "sans",
  label_fontface = "bold",
  hjust = -0.2,
  vjust = 1.1,
  axis = "l",
  rel_widths = c(1, 1, 1), 
  rel_heights = c(1, 1) 
) +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )
figure2
getwd()
ggsave("figure2_renew.png", figure2, width = 12, height = 12, dpi = 300)
ggsave("figure2.pdf", figure2, width = 18, height = 10)

