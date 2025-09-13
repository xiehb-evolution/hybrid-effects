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

###################################################################################################################
# FIGURE 1-AB: FST, Lambda, and Recombination Relationships
###################################################################################################################
# Panel A: FST vs Lambda relationship with recombination rate
sql = "SELECT a.fst_bin,(max(b.weighted_fst)+min(b.weighted_fst))/2 as fst,
       sum(b.recombinations)/578/2/count(*) as rec,
       sum(b.lambda)/count(*) as lambda
       FROM fst100k_autosomes a,renew_complete_data_with_overall b
       where a.chr=b.chr and a.window=b.window and b.Sex='Overall'
       group by a.fst_bin"
fst_lambda_data = dbGetQuery(con, sql)
# Calculate correlation for annotation
fst_lambda_cor = cor.test(fst_lambda_data$fst, fst_lambda_data$lambda)

panel_2a <- ggplot(fst_lambda_data, aes(x = fst, y = lambda)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#F8F9FA", color = NA),
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.grid.major = element_line(color = "#E9ECEF", size = 0.2),
    panel.grid.minor = element_line(color = "#E9ECEF", size = 0.1)
  ) +
  geom_smooth(method = "loess", se = TRUE, 
              color = "#495057", 
              fill = "#ADB5BD",
              linetype = "dashed", 
              size = 0.7, 
              span = 0.75) + 
  geom_hline(yintercept = mean(fst_lambda_data$lambda), 
             linetype = "dotted", color = "#6C757D", size = 0.4) +
  geom_vline(xintercept = mean(fst_lambda_data$fst), 
             linetype = "dotted", color = "#6C757D", size = 0.4) +
  geom_point(aes(size = rec, color = rec), alpha = 0.85) +
  scale_color_gradientn(
    colors = c("#F28C28", "#F7AD19", "#DF7861", "#A84448", "#772D53", "#541E54", "#371546", "#1A0B36", "#050225"),
    name = "Recombination rate",
    guide = "none"
  ) + 
  scale_size_continuous(range = c(2.5, 5), name = "Recombination rate") +
  labs(
    x = expression(F[ST]), 
    y = expression(lambda)
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 11, color = "#212529"),
    axis.text = element_text(size = 9, color = "#495057"),
    plot.title = element_text(face = "bold", size = 12, color = "#212529"),
    plot.subtitle = element_text(size = 9, color = "#6C757D", margin = margin(b = 15)),
    axis.line = element_line(color = "#343A40", size = 0.3)
  ) +
  sci_theme()+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  ) +
  guides(
    size = guide_legend(
      title.position = "top",
      override.aes = list(
        color = c("#F28C28", "#F7AD19", "#DF7861", "#A84448", "#541E54")
      )
    )
  )
panel_2a



# Panel B: Recombination rate vs Lambda
# Query data for recombination vs lambda analysis
sql = "SELECT a.fst_bin,(max(b.weighted_fst)+min(b.weighted_fst))/2 as fst,
       sum(b.recombinations)/578/2/count(*) as rec,
       sum(b.lambda)/count(*) as lambda
       FROM fst100k_autosomes a,renew_complete_data_with_overall b
       where a.chr=b.chr and a.window=b.window and b.Sex='Overall'
       group by a.fst_bin"
plot1data = dbGetQuery(con, sql)

# Calculate correlation for annotation
rec_lambda_cor <- cor.test(plot1data$rec, plot1data$lambda)
rec_lambda_lm <- lm(lambda ~ rec, data = plot1data)

# Create recombination vs lambda plot with linear model
panel_2b <- ggplot(plot1data, aes(x = rec, y = lambda)) +
  geom_vline(xintercept = mean(plot1data$rec, na.rm = TRUE), 
             color = "gray80", linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = mean(plot1data$lambda, na.rm = TRUE), 
             color = "gray80", linetype = "dashed", size = 0.3) +
  geom_smooth(method = "lm", color = "#E69F00", 
              alpha = 0.2, se = TRUE, size = 0.8) +
  geom_point(color = "#0072B2", alpha = 0.8, size = 3.5) +
  labs(
    x = "Recombination rate",
    y = expression(lambda),
    title = "Relationship between recombination rate and hybrid effect"
  ) +
  annotate("text", x = min(plot1data$rec, na.rm = TRUE), 
           y = max(plot1data$lambda, na.rm = TRUE),
           label = sprintf("r = %.3f, %s", 
                           rec_lambda_cor$estimate, 
                           format_pvalue(rec_lambda_cor$p.value)),
           hjust = 0, vjust = 1, size = 3.5) +
  sci_theme()
panel_2b


###########################################################################################################
# FIGURE 1-CD: Genomic Distribution and Characteristics of Hybrid Effects
###########################################################################################################
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

# Panel C: Recombination rate by hybrid effect type (log scale)
panel_3c <- ggplot(hybrid_effect_data, aes(x = tag, y = recombination_rate_all, fill = tag)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Recombination rate (log scale)",
    title = "Recombination rates across hybrid effect types"
  ) +
  scale_x_discrete(labels = c("Inbreeding\ndepression", "Hybrid\nvigor", "Hybrid\ndepression")) +
  sci_theme() +
  theme(legend.position = "none")

# Add statistical comparison
panel_3c <- panel_3c +
  stat_compare_means(comparisons = list(
    c("inbreeding depression", "hybrid vigor"),
    c("hybrid vigor", "hybrid depression"),
    c("inbreeding depression", "hybrid depression")),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    size = 3)
panel_3c



# Panel D: Genomic distribution of inbreeding depression
inbreeding_data <- hybrid_effect_data %>% 
  filter(tag == "inbreeding depression", chr >= 1, chr <= 18) %>%
  arrange(desc(lambda))
inbreeding_data$chr <- factor(inbreeding_data$chr, levels = 1:18)

panel_3d <- ggplot(inbreeding_data, aes(x = chr, y = window)) +
  geom_tile(aes(fill = lambda), width = 0.8) +
  scale_fill_viridis_c(option = "magma", 
                       name = expression(lambda), 
                       direction = 1) +
  labs(
    x = "Chromosome",
    y = "Window position",
    title = "Genomic distribution of inbreeding depression"
  ) +
  sci_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
panel_3d


# Panel E: Genomic distribution of hybrid depression
hybrid_dep_data <- hybrid_effect_data %>% 
  filter(tag == "hybrid depression", chr >= 1, chr <= 18) %>%
  arrange(desc(lambda))
hybrid_dep_data$chr <- factor(hybrid_dep_data$chr, levels = 1:18)

panel_3e <- ggplot(hybrid_dep_data, aes(x = chr, y = window)) +
  geom_tile(aes(fill = lambda), width = 0.8) +
  scale_fill_viridis_c(option = "magma", 
                       name = expression(lambda), 
                       direction = 1) +
  labs(
    x = "Chromosome",
    y = "Window position",
    title = "Genomic distribution of hybrid depression"
  ) +
  sci_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
panel_3e


panel_3d <- panel_3d + 
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(10, 10, 10, 10),
    legend.box.margin = margin(5, 5, 5, 5),
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.2),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )
panel_3e <- panel_3e + 
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(10, 10, 10, 10),
    legend.box.margin = margin(5, 5, 5, 5),
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.2),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

####################################################################################################################
panel_2a <- panel_2a + labs(title = NULL)
panel_2a
panel_2b <- panel_2b + labs(title = NULL)
panel_2b
panel_3d <- panel_3d + labs(title = NULL)
panel_3d
panel_3e <- panel_3e + labs(title = NULL)
panel_3e
figure3 <- plot_grid(
  panel_2a, panel_2b,
  panel_3d, panel_3e,
  ncol = 2, 
  align = 'v',
  labels = c("A", "B", "C", "D"),
  label_size = 16,
  label_fontfamily = "sans",
  label_fontface = "bold",
  hjust = -0.2,
  vjust = 1.1
)
figure3
ggsave("figure1_new.pdf", figure3, width = 12, height = 10, device = cairo_pdf)
getwd()
