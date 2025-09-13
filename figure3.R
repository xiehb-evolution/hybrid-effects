library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scales)
setwd("D:\\heterosis\\xie")
getwd()

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

# SQL Query to get data from database
sql_query_from_file <- paste0("select c.fst_bin,a.lambda,sum(abs(b.malecount_diff+b.femalecount_diff))/count(*) as countdiff,
sum(case b.malecount_diff+b.femalecount_diff>0 and d.malecount_diff+d.femalecount_diff>0 when 1 then 1 end) as count_increase,
sum(case b.malecount_diff+b.femalecount_diff<0 and d.malecount_diff+d.femalecount_diff<0 when 1 then 1 end) as count_decrease,
sum(case (b.malecount_diff+b.femalecount_diff)*(d.malecount_diff+d.femalecount_diff)<0 when 1 then 1 end) as count_conflicts,
sum(case b.malecount_diff+b.femalecount_diff=0 and d.malecount_diff+d.femalecount_diff=0 when 1 then 1 end) as count_nochange,
sum(case b.malecount_diff>0 and d.malecount_diff>0 when 1 then 1 end) as count_increase_male,
sum(case b.malecount_diff<0 and d.malecount_diff<0 when 1 then 1 end) as count_decrease_male,
sum(case b.malecount_diff=0 and d.malecount_diff=0 when 1 then 1 end) as count_nochange_male,
sum(case b.malecount_diff>0 and d.malecount_diff>0 when 1 then 1 end)/sum(case b.malecount_diff<0 and d.malecount_diff<0 when 1 then 1 end) as maleratio,
sum(case b.femalecount_diff>0 and d.femalecount_diff>0 when 1 then 1 end) as count_increase_female,
sum(case b.femalecount_diff<0 and d.femalecount_diff<0 when 1 then 1 end) as count_decrease_female,
sum(case b.femalecount_diff=0 and d.femalecount_diff=0 when 1 then 1 end) as count_nochange_female,
sum(case b.femalecount_diff>0 and d.femalecount_diff>0 when 1 then 1 end)/sum(case b.femalecount_diff<0 and d.femalecount_diff<0 when 1 then 1 end) as femaleratio
from renew_complete_data_with_overall a,window100k_single_site_stat_mutant b,fst100k_autosomes c,window100k_single_site_stat_mutant d
where a.sex='Overall' and a.chr=c.chr and a.window=c.window 
and a.chr=b.chr and a.window=b.window and b.paternalinheritance=b.maternalinheritance and b.paternalinheritance=0
and a.chr=d.chr and a.window=d.window and d.paternalinheritance=d.maternalinheritance and d.paternalinheritance=1 
and b.maternalinheritance_mutant=d.maternalinheritance_mutant and b.paternalinheritance_mutant=d.paternalinheritance_mutant
group by c.fst_bin,a.lambda
having count_increase>=10
limit 100")
result <- dbGetQuery(con, sql_query_from_file)
# Get lambda data for males and females
male_query <- "
  SELECT fst_bin, AVG(lambda) as male_lambda
  FROM renew_complete_data_with_overall
  WHERE sex = 'Male'
  GROUP BY fst_bin
"
male_lambda <- dbGetQuery(con, male_query)
female_query <- "
  SELECT fst_bin, AVG(lambda) as female_lambda
  FROM renew_complete_data_with_overall
  WHERE sex = 'Female'
  GROUP BY fst_bin
"
female_lambda <- dbGetQuery(con, female_query)
# Merge data together
result_enhanced <- merge(result, male_lambda, by = "fst_bin", all.x = TRUE)
result_enhanced <- merge(result_enhanced, female_lambda, by = "fst_bin", all.x = TRUE)
# Get FST range information
fst_start_query <- "
  SELECT fst_bin, AVG(fst_start) as fst_start, AVG(fst_end) as fst_end
  FROM renew_complete_data_with_overall
  WHERE sex = 'Overall'
  GROUP BY fst_bin
"
fst_start_result <- dbGetQuery(con, fst_start_query)
result_enhanced <- merge(result_enhanced, fst_start_result, by = "fst_bin", all.x = TRUE)
# Calculate additional metrics
result_enhanced <- result_enhanced %>%
  mutate(
    het_advantage_ratio = count_increase / count_decrease,
    fst_midpoint = (fst_start + fst_end) / 2,
    het_advantage_ratio_male = count_increase_male / count_decrease_male,
    het_advantage_ratio_female = count_increase_female / count_decrease_female
  )

dbExecute(con, "
CREATE TABLE IF NOT EXISTS hybrid_effect_analysis (
    fst_bin INT,
    lambda FLOAT,
    countdiff FLOAT,
    count_increase INT,
    count_decrease INT,
    count_conflicts INT,
    count_nochange INT,
    count_increase_male INT,
    count_decrease_male INT,
    count_nochange_male INT,
    maleratio FLOAT,
    count_increase_female INT,
    count_decrease_female INT,
    count_nochange_female INT,
    femaleratio FLOAT,
    male_lambda FLOAT,
    female_lambda FLOAT,
    fst_start FLOAT,
    fst_end FLOAT,
    het_advantage_ratio FLOAT,
    fst_midpoint FLOAT,
    het_advantage_ratio_male FLOAT,
    het_advantage_ratio_female FLOAT,
    PRIMARY KEY (fst_bin, lambda)
);")
#dbWriteTable(con, "hybrid_effect_analysis", result_enhanced, append = TRUE, row.names = FALSE)

#===============================#
# Figure 3B: Hump-shaped distribution
#===============================#

male_data <- data.frame(
  fst_bin = result_enhanced$fst_bin,
  fst_midpoint = result_enhanced$fst_midpoint,
  ratio = result_enhanced$maleratio,
  sex = "Male"
)
female_data <- data.frame(
  fst_bin = result_enhanced$fst_bin,
  fst_midpoint = result_enhanced$fst_midpoint,
  ratio = result_enhanced$femaleratio,
  sex = "Female"
)
result_long <- rbind(male_data, female_data)
publication_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Arial", color = "#2c3e50"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    panel.grid.major = element_line(color = "#ecf0f1", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "#2c3e50", size = 0.3),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "#bdc3c7", size = 0.3),
    legend.key = element_rect(fill = NA),
    legend.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 9)
  )
result_enhanced$Group <- "All"
result_male <- result_long %>% 
  filter(sex == "Male") %>%
  rename(het_advantage_ratio = ratio) %>%
  mutate(Group = "Male")

result_female <- result_long %>% 
  filter(sex == "Female") %>%
  rename(het_advantage_ratio = ratio) %>%
  mutate(Group = "Female")

combined_data <- bind_rows(
  result_enhanced,
  result_male,
  result_female
)

combined_data$Group <- factor(combined_data$Group, levels = c("All", "Male", "Female"))

plot3B <- ggplot(combined_data, aes(x = fst_midpoint, y = het_advantage_ratio, color = Group)) +
  geom_point(aes(shape = Group), alpha = 0.8, size = 2) +
  geom_smooth(aes(fill = Group), method = "loess", se = TRUE, 
              alpha = 0.15, size = 0.8, span = 0.95) +
  geom_hline(yintercept = 1, linetype = "dotted", 
             color = "#7f8c8d", size = 0.5) +
  scale_color_manual(values = c(
    "All" = "#27ae60",    # 绿色
    "Male" = "#3366CC",   # 深蓝色
    "Female" = "#CC3366"  # 粉红色
  )) +
  scale_fill_manual(values = c(
    "All" = "#27ae60",    # 绿色
    "Male" = "#3366CC",   # 深蓝色
    "Female" = "#CC3366"  # 粉红色
  )) +
  scale_shape_manual(values = c(
    "All" = 0,      # 空心圆形
    "Male" = 19,    # 实心圆形
    "Female" = 17   # 实心三角形
  )) +
  scale_y_continuous(
    trans = "log2", 
    breaks = c(0.25, 0.5, 1, 2, 4, 8), 
    labels = c("1/4", "1/2", "1", "2", "4", "8")
  ) +
  labs(
    x = expression(F[ST]), 
    y = "Heterozygote advantage/disadvantage ratio"
  ) +
  clean_theme+sci_theme(11)+
  # Position the legend inside the top-right corner of the plot
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.2)
  )
plot3B


#===============================#
# Figure 3C
#===============================#
# Calculate correlation coefficient
cor_result <- cor.test(result_enhanced$lambda, result_enhanced$het_advantage_ratio, 
                       method = "pearson", use = "complete.obs")
cor_value <- cor_result$estimate
p_value <- cor_result$p.value

cor_text <- sprintf("r = %.3f, p = %.3e", cor_value, p_value)

plot3C <- ggplot(result_enhanced, 
                          aes(x = lambda, y = het_advantage_ratio)) +

  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", size = 0.5) +
  
  geom_point(color = "#2980b9", size = 2, alpha = 0.7) +
  
  geom_smooth(method = "lm", color = "#c0392b", 
              fill = "#c0392b", alpha = 0.1, size = 0.7) +
  
  scale_y_continuous(
    trans = "log2", 
    breaks = c(0.25, 0.5, 1, 2, 4, 8), 
    labels = c("1/4", "1/2", "1", "2", "4", "8")
  ) +
  
  labs(
    x = expression(lambda),
    y = "Heterozygote advantage/disadvantage ratio"
  ) +
  
  annotate("text", 
           x = min(result_enhanced$lambda, na.rm = TRUE), 
           y = max(result_enhanced$het_advantage_ratio, na.rm = TRUE),
           label = cor_text,
           hjust = 0, vjust = 1, size = 3.5) +
  
  theme_minimal() +sci_theme(11)+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "whitesmoke"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )
plot3C



library(cowplot)
figure3 <- plot_grid(
  NULL, plot3B,plot3C,
  ncol = 3,
  align = 'hv',
  labels = c("A", "B", "C"),  # B位置留空标签
  label_size = 16,
  label_fontfamily = "sans",
  label_fontface = "bold",
  hjust = -0.2,
  vjust = 1.1,
  axis = "l",
  rel_widths = c(1, 1, 1),  # 三列宽度相同
  rel_heights = c(1, 1)  # 两行高度相同
) +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )
figure3
ggsave("figure3.pdf", figure3, width = 18, height = 5)
