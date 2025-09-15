# ====================================================
# Build Script for Integrated Genetic Data Analysis
# This script creates and populates a comprehensive table
# integrating genetic diversity, recombination, and FST data
# ====================================================

# Set working directory and load required libraries
setwd("E:\\xiehb_sync\\F2\\data\\sscrofa11.data")

require(RMySQL)
require(DBI)
require(dplyr)

# Database connection parameters
mysqlserver <- "10.0.16.70"
mysqlport <- 3306
mysqldbname <- "speciation"
mysqluser <- "xiehb"
mysqlpassword <- "password"

# Establish database connection
con <- dbConnect(dbDriver("MySQL"), 
                host = mysqlserver, 
                port = mysqlport, 
                dbname = mysqldbname, 
                user = mysqluser, 
                password = mysqlpassword)

#--------------------------------------------------------------------
# Create the base table structure for integrated genetic analysis
#--------------------------------------------------------------------
create_table_sql <- "
CREATE TABLE IF NOT EXISTS renew_complete_data_with_overall (
    chr INT NOT NULL COMMENT 'Chromosome number',
    window INT NOT NULL COMMENT 'Genomic window position',
    lambda DECIMAL(10,4) COMMENT 'Genetic diversity measure',
    WEIGHTED_FST DECIMAL(10,8) COMMENT 'Population differentiation index',
    snps INT COMMENT 'Number of SNPs in the window',
    recombination_rate DECIMAL(20,18) COMMENT 'Recombination rate in the window',
    Sex ENUM('Female', 'Male', 'Overall') COMMENT 'Sex classification',
    PRIMARY KEY (chr, window, Sex)
) COMMENT 'Integrated table with genetic diversity, recombination, and FST data'"

dbSendQuery(con, create_table_sql)

#--------------------------------------------------------------------
# Populate the table with base data from various sources
#--------------------------------------------------------------------
populate_base_data_sql <- "
INSERT INTO renew_complete_data_with_overall (chr, window, lambda, WEIGHTED_FST, snps, recombination_rate, Sex)
SELECT 
    c.chr,
    c.window,
    c.lambda,
    c.WEIGHTED_FST,
    c.snps,
    w.recombination_rate,
    c.Sex
FROM 
    complete_data_with_overall_table c
JOIN 
    window100k_recombination_statistics w ON c.chr = w.chr AND c.window = w.window"

dbSendQuery(con, populate_base_data_sql)

#--------------------------------------------------------------------
# Add recombined_count column and populate it
#--------------------------------------------------------------------
# Add the column
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall
                 ADD COLUMN recombined_count INT COMMENT 'Count of recombination events in the window'")

# Populate recombined_count values
update_recombined_count_sql <- "
UPDATE renew_complete_data_with_overall AS c
JOIN window100k_recombination_statistics AS w
  ON c.chr = w.chr AND c.window = w.window
SET c.recombined_count = w.recombined_count"

dbSendQuery(con, update_recombined_count_sql)

#--------------------------------------------------------------------
# Update recombination_rate with higher precision values
#--------------------------------------------------------------------
# Modify column precision
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall 
                 MODIFY COLUMN recombination_rate DECIMAL(20,18)")

# Update with high precision recombination rates
update_recombination_rate_sql <- "
UPDATE renew_complete_data_with_overall AS c
JOIN window100k_recombination_statistics AS w
  ON c.chr = w.chr AND c.window = w.window
SET c.recombination_rate = w.recombination_rate"

dbSendQuery(con, update_recombination_rate_sql)

#--------------------------------------------------------------------
# Add and populate mean_rec_rate (average recombination rate by sex and lambda)
#--------------------------------------------------------------------
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall 
                 ADD COLUMN mean_rec_rate DECIMAL(20,18) COMMENT 'Average recombination rate by sex and lambda'")

calculate_mean_rec_rate_sql <- "
UPDATE renew_complete_data_with_overall AS t1
JOIN (
    SELECT 
        Sex, 
        lambda, 
        AVG(recombination_rate) AS avg_rate
    FROM 
        renew_complete_data_with_overall
    GROUP BY 
        Sex, lambda
) AS t2
ON t1.Sex = t2.Sex AND t1.lambda = t2.lambda
SET t1.mean_rec_rate = t2.avg_rate"

dbSendQuery(con, calculate_mean_rec_rate_sql)

#--------------------------------------------------------------------
# Add and populate tag column based on FST quantiles
#--------------------------------------------------------------------
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall 
                 ADD COLUMN tag VARCHAR(50) DEFAULT NULL COMMENT 'Classification based on FST quantiles'")

# Classification based on FST quantiles:
# 0% ~ 25% quantile: inbreeding depression
# 25% ~ 75% quantile: hybrid vigor  
# 75% ~ 100% quantile: hybrid depression
classify_by_fst_sql <- "
UPDATE renew_complete_data_with_overall
SET tag = 
    CASE 
        WHEN WEIGHTED_FST < 0.09454525 THEN 'inbreeding depression' 
        WHEN WEIGHTED_FST BETWEEN 0.09454525 AND 0.23909900 THEN 'hybrid vigor'
        WHEN WEIGHTED_FST > 0.23909900 THEN 'hybrid depression'
        ELSE NULL
    END"

dbSendQuery(con, classify_by_fst_sql)

# Alternative classification method based on lambda and FST values
alternative_classify_sql <- "
UPDATE renew_complete_data_with_overall
SET tag = 
    CASE 
        WHEN lambda < 33.5 AND WEIGHTED_FST < 0.2 THEN 'inbreeding depression'
        WHEN lambda < 33.5 AND WEIGHTED_FST > 0.2 THEN 'hybrid depression'
        WHEN lambda > 33.5 THEN 'hybrid vigor'
        ELSE NULL
    END"

# Uncomment the line below to use alternative classification
# dbSendQuery(con, alternative_classify_sql)

#--------------------------------------------------------------------
# Add and populate sort_rec_rate (average recombination rate by sex and tag)
#--------------------------------------------------------------------
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall 
                 ADD COLUMN sort_rec_rate DECIMAL(20,18) COMMENT 'Average recombination rate by sex and tag category'")

calculate_sort_rec_rate_sql <- "
UPDATE renew_complete_data_with_overall t1
JOIN (
    SELECT 
        Sex, 
        tag, 
        AVG(recombination_rate) AS avg_rec_rate
    FROM 
        renew_complete_data_with_overall
    GROUP BY 
        Sex, tag
) t2 ON t1.Sex = t2.Sex AND t1.tag = t2.tag
SET t1.sort_rec_rate = t2.avg_rec_rate"

dbSendQuery(con, calculate_sort_rec_rate_sql)

#--------------------------------------------------------------------
# Add FST start and end positions
#--------------------------------------------------------------------
dbSendQuery(con, "ALTER TABLE renew_complete_data_with_overall
                 ADD COLUMN fst_start DECIMAL(10,8) COMMENT 'Start position of FST window',
                 ADD COLUMN fst_end DECIMAL(10,8) COMMENT 'End position of FST window'")

update_fst_positions_sql <- "
UPDATE renew_complete_data_with_overall AS t1
JOIN all_data AS t2
ON t1.chr = t2.chr AND t1.window = t2.window
SET t1.fst_start = t2.fst_start,
    t1.fst_end = t2.fst_end"

dbSendQuery(con, update_fst_positions_sql)

#--------------------------------------------------------------------
# Create indexes for better query performance
#--------------------------------------------------------------------
index_queries <- c(
  "CREATE INDEX idx_renew_sex_tag ON renew_complete_data_with_overall(Sex, tag)",
  "CREATE INDEX idx_renew_fst ON renew_complete_data_with_overall(WEIGHTED_FST)",
  "CREATE INDEX idx_renew_lambda ON renew_complete_data_with_overall(lambda)",
  "CREATE INDEX idx_renew_chr_window ON renew_complete_data_with_overall(chr, window)"
)

for(index_sql in index_queries) {
  tryCatch({
    dbSendQuery(con, index_sql)
  }, error = function(e) {
    # Index might already exist, continue
  })
}

#--------------------------------------------------------------------
# Export final integrated data table for further analysis
#--------------------------------------------------------------------
# Query the completed table
final_data <- dbGetQuery(con, "SELECT * FROM renew_complete_data_with_overall ORDER BY chr, window, Sex")

# Save to local file for backup and further R analysis
write.table(final_data, "integrated_genetic_data_complete.txt", 
           quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Close database connection
dbDisconnect(con)

# Display summary statistics
summary_stats <- final_data %>%
  group_by(Sex, tag) %>%
  summarise(
    count = n(),
    mean_lambda = mean(lambda, na.rm = TRUE),
    mean_fst = mean(WEIGHTED_FST, na.rm = TRUE),
    mean_recombination = mean(recombination_rate, na.rm = TRUE),
    .groups = 'drop'
  )

print("Summary statistics by Sex and Tag:")
print(summary_stats)

#--------------------------------------------------------------------
# Analysis Tables for Heterosis Studies
#--------------------------------------------------------------------
# Create specialized analysis tables for studying bidirectional heterosis
# and deviation patterns within FST bins

# Create table to analyze bidirectional heterosis (BPH) statistics within FST bins
create_bph_table_sql <- "
CREATE TABLE stat_bph_in_fst_bin 
SELECT c.fst_bin, a.lambda,
    -- Count bidirectional heterosis cases (both male and female trait differences are non-zero and same direction)
    SUM(CASE WHEN b.maletrait_diff*b.maletrait_diff>0 AND d.maletrait_diff*d.femaletrait_diff>0 THEN 1 ELSE 0 END) as count_bph,
    -- Count non-BPH cases
    SUM(CASE WHEN NOT (b.maletrait_diff*b.maletrait_diff>0 AND d.maletrait_diff*d.femaletrait_diff>0) THEN 1 ELSE 0 END) as count_non_bph,
    -- Count male-specific BPH cases (male trait differences have same direction)
    SUM(CASE WHEN b.maletrait_diff*d.maletrait_diff>0 THEN 1 ELSE 0 END) as count_bph_male,
    -- Count male-specific non-BPH cases (male trait differences have opposite directions)
    SUM(CASE WHEN b.maletrait_diff*d.maletrait_diff<0 THEN 1 ELSE 0 END) as count_non_bph_male,
    -- Count female-specific BPH cases (female trait differences have same direction)
    SUM(CASE WHEN b.femaletrait_diff*d.femaletrait_diff>0 THEN 1 ELSE 0 END) as count_bph_female,
    -- Count female-specific non-BPH cases (female trait differences have opposite directions)
    SUM(CASE WHEN b.femaletrait_diff*d.femaletrait_diff<0 THEN 1 ELSE 0 END) as count_non_bph_female
FROM renew_complete_data_with_overall a,
     window100k_single_site_trait_stat_mutant_deviation_from_mean b,
     fst100k_autosomes c,
     window100k_single_site_trait_stat_mutant_deviation_from_mean d
WHERE a.sex='Overall' AND a.chr=c.chr AND a.window=c.window 
    AND a.chr=b.chr AND a.window=b.window AND b.paternalinheritance=b.maternalinheritance AND b.paternalinheritance=0
    AND a.chr=d.chr AND a.window=d.window AND d.paternalinheritance=d.maternalinheritance AND d.paternalinheritance=1 
    AND b.maternalinheritance_mutant=d.maternalinheritance_mutant AND b.paternalinheritance_mutant=d.paternalinheritance_mutant
GROUP BY c.fst_bin, a.lambda
LIMIT 100"

tryCatch({
  dbSendQuery(con, "DROP TABLE IF EXISTS stat_bph_in_fst_bin")
  dbSendQuery(con, create_bph_table_sql)
}, error = function(e) {
  print(paste("Warning: Could not create BPH analysis table:", e$message))
})

# Create table to calculate Mid-Parent Heterosis (MPH) within FST bins
create_mph_table_sql <- "
CREATE TABLE stat_mph_in_fst_bin 
SELECT c.fst_bin,
    -- Calculate average MPH for male traits
    SUM((2*b.maletrait_diff + b.maletrait - d.maletrait)/(b.maletrait+d.maletrait))/COUNT(*) as mph_male,
    -- Calculate average MPH for female traits
    SUM((b.femaletrait_diff*2 + b.femaletrait - d.femaletrait)/(b.femaletrait+d.femaletrait))/COUNT(*) as mph_female
FROM renew_complete_data_with_overall a,
     window100k_single_site_trait_stat_mutant_deviation_from_mean b,
     fst100k_autosomes c,
     window100k_single_site_trait_stat_mutant_deviation_from_mean d
WHERE a.sex='Overall' AND a.chr=c.chr AND a.window=c.window 
    AND a.chr=b.chr AND a.window=b.window AND b.paternalinheritance=b.maternalinheritance AND b.paternalinheritance=0
    AND a.chr=d.chr AND a.window=d.window AND d.paternalinheritance=d.maternalinheritance AND d.paternalinheritance=1 
    AND b.maternalinheritance_mutant=d.maternalinheritance_mutant AND b.paternalinheritance_mutant=d.paternalinheritance_mutant
GROUP BY c.fst_bin
LIMIT 100"

tryCatch({
  dbSendQuery(con, "DROP TABLE IF EXISTS stat_mph_in_fst_bin")
  dbSendQuery(con, create_mph_table_sql)
}, error = function(e) {
  print(paste("Warning: Could not create MPH analysis table:", e$message))
})

# Create table to analyze deviation from sex-specific means within FST bins (version 1)
create_dev_table2_sql <- "
CREATE TABLE stat_dev_from_sex_mean_in_fst_bin2
SELECT c.fst_bin, a.lambda, COUNT(*) as count,
    -- Calculate relative difference in male trait values between genotypes
    SUM(ABS(b.maletrait+2*b.maletrait_diff-d.maletrait)/(b.maletrait+d.maletrait))/COUNT(*) as malediff,
    -- Calculate relative difference in female trait values between genotypes
    SUM(ABS(b.femaletrait+2*b.femaletrait_diff-d.femaletrait)/(b.femaletrait+d.femaletrait))/COUNT(*) as femalediff
FROM renew_complete_data_with_overall a,
     window100k_single_site_trait_stat_mutant_deviation_from_mean b,
     fst100k_autosomes c,
     window100k_single_site_trait_stat_mutant_deviation_from_mean d
WHERE a.sex='Overall' AND a.chr=c.chr AND a.window=c.window 
    AND a.chr=b.chr AND a.window=b.window AND b.paternalinheritance=b.maternalinheritance AND b.paternalinheritance=0
    AND a.chr=d.chr AND a.window=d.window AND d.paternalinheritance=d.maternalinheritance AND d.paternalinheritance=1 
    AND b.maternalinheritance_mutant=d.maternalinheritance_mutant AND b.paternalinheritance_mutant=d.paternalinheritance_mutant
    AND b.trait_id=d.trait_id
GROUP BY c.fst_bin, a.lambda
HAVING COUNT(*)>1000
LIMIT 100"

tryCatch({
  dbSendQuery(con, "DROP TABLE IF EXISTS stat_dev_from_sex_mean_in_fst_bin2")
  dbSendQuery(con, create_dev_table2_sql)
}, error = function(e) {
  print(paste("Warning: Could not create deviation analysis table (version 2):", e$message))
})

# Create table to analyze deviation from sex-specific means within FST bins (version 2 - improved)
create_dev_table_sql <- "
CREATE TABLE stat_dev_from_sex_mean_in_fst_bin 
SELECT c.fst_bin, a.lambda, COUNT(*) as count,
    -- Calculate absolute difference in male deviation values between genotypes
    SUM(ABS(b.maledev-d.maledev))/COUNT(*) as malediff,
    -- Calculate absolute difference in female deviation values between genotypes
    SUM(ABS(b.femaledev-d.femaledev))/COUNT(*) as femalediff
FROM renew_complete_data_with_overall a,
     window100k_single_site_trait_stat_mutant_deviation_from_mean b,
     fst100k_autosomes c,
     window100k_single_site_trait_stat_mutant_deviation_from_mean d
WHERE a.sex='Overall' AND a.chr=c.chr AND a.window=c.window 
    AND a.chr=b.chr AND a.window=b.window AND b.paternalinheritance=b.maternalinheritance AND b.paternalinheritance=0
    AND a.chr=d.chr AND a.window=d.window AND d.paternalinheritance=d.maternalinheritance AND d.paternalinheritance=1 
    AND b.maternalinheritance_mutant=d.maternalinheritance_mutant AND b.paternalinheritance_mutant=d.paternalinheritance_mutant
    AND b.trait_id=d.trait_id
GROUP BY c.fst_bin, a.lambda
HAVING COUNT(*)>1000
LIMIT 100"

tryCatch({
  dbSendQuery(con, "DROP TABLE IF EXISTS stat_dev_from_sex_mean_in_fst_bin")
  dbSendQuery(con, create_dev_table_sql)
}, error = function(e) {
  print(paste("Warning: Could not create deviation analysis table:", e$message))
})

#--------------------------------------------------------------------
# Export heterosis analysis tables for further R analysis
#--------------------------------------------------------------------
# Export BPH analysis results if table exists
tryCatch({
  bph_data <- dbGetQuery(con, "SELECT * FROM stat_bph_in_fst_bin ORDER BY fst_bin, lambda")
  write.table(bph_data, "bidirectional_heterosis_by_fst_bin.txt", 
             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  print(paste("Exported", nrow(bph_data), "BPH analysis records"))
}, error = function(e) {
  print("BPH analysis table not available for export")
})

# Export MPH analysis results if table exists
tryCatch({
  mph_data <- dbGetQuery(con, "SELECT * FROM stat_mph_in_fst_bin ORDER BY fst_bin")
  write.table(mph_data, "mid_parent_heterosis_by_fst_bin.txt", 
             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  print(paste("Exported", nrow(mph_data), "MPH analysis records"))
}, error = function(e) {
  print("MPH analysis table not available for export")
})

# Export deviation analysis results if table exists
tryCatch({
  dev_data <- dbGetQuery(con, "SELECT * FROM stat_dev_from_sex_mean_in_fst_bin ORDER BY fst_bin, lambda")
  write.table(dev_data, "sex_deviation_analysis_by_fst_bin.txt", 
             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  print(paste("Exported", nrow(dev_data), "deviation analysis records"))
}, error = function(e) {
  print("Deviation analysis table not available for export")
})

# ====================================================
# End of Build Script - Integrated Genetic Data Table
# ====================================================
