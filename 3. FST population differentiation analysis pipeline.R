#------------------------------------------------------------------------------------------------------------
# FST Population Differentiation Analysis and Data Import Pipeline
# This script calculates FST (Fixation Index) values using VCFtools and imports results into MySQL database
# FST measures population differentiation due to genetic structure and is commonly used in evolutionary studies
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
# FST Calculation Using VCFtools (Command Line Section)
# The following commands should be executed in the terminal before running this R script
# These commands calculate FST values in 100kb sliding windows across the genome
#------------------------------------------------------------------------------------------------------------

# Load VCFtools software environment
# source /public/software/profile.d/vcftools.sh

# Navigate to working directory containing VCF files and population information
# cd /public/home/xiehaibing1/12.snp/f0

# Calculate FST values in 100kb windows between two populations
# --vcf: Input VCF file containing genotype data for all individuals
# --weir-fst-pop: Population files containing sample IDs for each population
# --fst-window-size: Window size in base pairs (100000 = 100kb)
# --fst-window-step: Step size for sliding windows (50000 = 50kb overlap)
# --out: Output file prefix for results

# vcftools --vcf f0.vcf \
#   --weir-fst-pop population1.txt \
#   --weir-fst-pop population2.txt \
#   --fst-window-size 100000 \
#   --fst-window-step 50000 \
#   --out min_lw_fst100k_results

# Alternative: Calculate per-SNP FST values (for comparison or detailed analysis)
# vcftools --vcf f0.vcf \
#   --weir-fst-pop population1.txt \
#   --weir-fst-pop population2.txt \
#   --out fst_snp_results

# Expected output files:
# - min_lw_fst100k_results.windowed.weir.fst (windowed FST results)
# - min_lw_fst100k_results.log (analysis log file)

#------------------------------------------------------------------------------------------------------------

# Load required libraries for database connectivity and data processing
library(RMySQL)  # MySQL database interface for R
library(readr)   # Fast and friendly data reading functions

#------------------------------------------------------------------------------------------------------------
# File Input and Data Loading Section
# Set working directory and load FST calculation results from windowed analysis
#------------------------------------------------------------------------------------------------------------

# Set working directory containing FST analysis results
setwd("F:\\data")

# Read FST data from VCFtools windowed Weir's FST output file
# This file contains FST values calculated in 100kb sliding windows across the genome
# File format: tab-delimited with headers (CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST, MEAN_FST)
fst_data <- read_delim("min_lw_fst100k_results.windowed.weir.fst", 
                       delim = "\t",           # Tab-separated values
                       col_names = TRUE,       # First row contains column names
                       trim_ws = TRUE)         # Remove leading/trailing whitespace

# Display first few rows to verify data structure
head(fst_data)

#------------------------------------------------------------------------------------------------------------
# Database Schema Preparation
# Analyze data characteristics to optimize database table structure
#------------------------------------------------------------------------------------------------------------

# Determine maximum chromosome identifier length for optimal VARCHAR sizing
# This ensures the database schema can accommodate all chromosome names without truncation
max_chrom_length <- max(nchar(as.character(fst_data$CHROM)))
print(paste("Maximum CHROM length:", max_chrom_length))

#------------------------------------------------------------------------------------------------------------
# Database Table Creation
# Create optimized table structure for FST data storage and analysis
#------------------------------------------------------------------------------------------------------------

# Remove existing table if present to ensure clean import
dbSendQuery(con, "DROP TABLE IF EXISTS fst100k_results")

# Create new table with appropriate data types and constraints
# Table structure optimized for genomic coordinate queries and FST value analysis
dbSendQuery(con, "CREATE TABLE fst100k_results (
              id INT AUTO_INCREMENT PRIMARY KEY,           -- Unique identifier for each record
              chrom VARCHAR(20),                           -- Chromosome identifier (e.g., '1', '2', 'X', 'Y')
              bin_start INT,                               -- Start position of 100kb window (bp)
              bin_end INT,                                 -- End position of 100kb window (bp)
              n_variants INT,                              -- Number of SNPs within the window
              weighted_fst DOUBLE,                         -- Weir's weighted FST estimate
              mean_fst DOUBLE                              -- Mean FST value across variants in window
            )")

#------------------------------------------------------------------------------------------------------------
# Data Import Process
# Insert FST data row by row with proper handling of numerical precision
#------------------------------------------------------------------------------------------------------------

# Iterate through each row of FST data for database insertion
# Using row-by-row insertion to handle potential data type conversion issues
for (i in 1:nrow(fst_data)) {
  # Construct SQL INSERT statement with proper value formatting
  # Note: String values are quoted, numerical values are unquoted
  query <- paste0("INSERT INTO fst100k_results (chrom, bin_start, bin_end, n_variants, weighted_fst, mean_fst) VALUES ('",
                  fst_data$CHROM[i], "', ",              # Chromosome as string
                  fst_data$BIN_START[i], ", ",           # Window start position
                  fst_data$BIN_END[i], ", ",             # Window end position
                  fst_data$N_VARIANTS[i], ", ",          # Variant count
                  fst_data$WEIGHTED_FST[i], ", ",        # Weighted FST value
                  fst_data$MEAN_FST[i], ")")             # Mean FST value
  
  # Execute the INSERT query
  dbSendQuery(con, query)
}

#------------------------------------------------------------------------------------------------------------
# Data Verification and Quality Control
# Verify successful data import and display sample records
#------------------------------------------------------------------------------------------------------------

# Query first 10 records to verify successful data import
result <- dbGetQuery(con, "SELECT * FROM fst100k_results LIMIT 10")
print(result)

#------------------------------------------------------------------------------------------------------------
# Database Optimization and Quality Control
# Create indexes for efficient querying and perform data validation
#------------------------------------------------------------------------------------------------------------

# Create indexes for optimized query performance
# Index on genomic coordinates for position-based queries
dbSendQuery(con, "CREATE INDEX idx_chrom_pos ON fst100k_results(chrom, bin_start, bin_end)")

# Index on FST values for value-based filtering and sorting
dbSendQuery(con, "CREATE INDEX idx_fst_values ON fst100k_results(weighted_fst, mean_fst)")

print("Database indexes created successfully")

#------------------------------------------------------------------------------------------------------------
# Data Quality Assessment
# Perform comprehensive quality checks on imported FST data
#------------------------------------------------------------------------------------------------------------

# Check for negative FST values (may indicate population structure issues or technical problems)
negative_fst_query <- "SELECT COUNT(*) as negative_count FROM fst100k_results WHERE weighted_fst < 0 OR mean_fst < 0"
negative_fst_result <- dbGetQuery(con, negative_fst_query)
print(paste("Number of windows with negative FST values:", negative_fst_result$negative_count))

# Verify genomic coordinate consistency (check for proper window boundaries)
coordinate_check_query <- "SELECT COUNT(*) as inconsistent_windows FROM fst100k_results WHERE bin_end <= bin_start"
coordinate_result <- dbGetQuery(con, coordinate_check_query)
print(paste("Number of windows with inconsistent coordinates:", coordinate_result$inconsistent_windows))

# Analyze FST distribution across chromosomes
chromosome_stats_query <- "SELECT chrom, 
                           COUNT(*) as window_count,
                           AVG(weighted_fst) as mean_weighted_fst,
                           MIN(weighted_fst) as min_fst,
                           MAX(weighted_fst) as max_fst,
                           AVG(n_variants) as mean_variants_per_window
                           FROM fst100k_results 
                           GROUP BY chrom 
                           ORDER BY chrom"
chromosome_stats <- dbGetQuery(con, chromosome_stats_query)
print("FST statistics by chromosome:")
print(chromosome_stats)

# Summary statistics for overall dataset
overall_stats_query <- "SELECT COUNT(*) as total_windows,
                        AVG(weighted_fst) as overall_mean_fst,
                        STDDEV(weighted_fst) as fst_stddev,
                        AVG(n_variants) as mean_variants,
                        MIN(n_variants) as min_variants,
                        MAX(n_variants) as max_variants
                        FROM fst100k_results"
overall_stats <- dbGetQuery(con, overall_stats_query)
print("Overall dataset statistics:")
print(overall_stats)

#------------------------------------------------------------------------------------------------------------
# Data Integration Preparation
# Prepare queries for integration with other genetic datasets
#------------------------------------------------------------------------------------------------------------

# Example query template for joining with QTL mapping results (commented for reference)
# qtl_fst_join_query <- "SELECT q.*, f.weighted_fst, f.mean_fst, f.n_variants
#                        FROM qtl_results q
#                        JOIN fst100k_results f ON q.chr = f.chrom 
#                        AND q.position >= f.bin_start 
#                        AND q.position <= f.bin_end
#                        WHERE f.weighted_fst > 0.1"  # Filter for high differentiation regions

# Example query for identifying high FST regions (potential selection signatures)
high_fst_regions_query <- "SELECT chrom, bin_start, bin_end, weighted_fst, mean_fst, n_variants
                          FROM fst100k_results 
                          WHERE weighted_fst > (SELECT AVG(weighted_fst) + 2*STDDEV(weighted_fst) FROM fst100k_results)
                          ORDER BY weighted_fst DESC"
high_fst_regions <- dbGetQuery(con, high_fst_regions_query)
print("High FST regions (>2 standard deviations above mean):")
print(head(high_fst_regions, 10))

#------------------------------------------------------------------------------------------------------------
# Analysis Complete
# FST data has been successfully calculated, imported, indexed, and validated
#------------------------------------------------------------------------------------------------------------
print("FST analysis pipeline completed successfully!")
print("Database table 'fst100k_results' is ready for downstream analysis.")

#------------------------------------------------------------------------------------------------------------