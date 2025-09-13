#--------------------------------------------------------------------
# 4.create_recombination_tables.R
# This script creates recombination statistics tables for genomic analysis
# It processes F2 inheritance data to identify recombination events
# Author: Modified from SQL version to R implementation
#--------------------------------------------------------------------

# Load required libraries
require(RMySQL)

# Database connection parameters (modify as needed)
mysqlserver = "XXXXX"
mysqlport = XXXXX
mysqldbname = "XXXXX"
mysqluser = "XXXXX"
mysqlpassword = "XXXXX"

# Establish database connection
con <- dbConnect(dbDriver("MySQL"), 
                 host = mysqlserver, 
                 port = mysqlport, 
                 dbname = mysqldbname, 
                 user = mysqluser, 
                 password = mysqlpassword)

#--------------------------------------------------------------------
# Create table to store recombination statistics for 100kb genomic windows
# This table aggregates recombination counts from allele count data, grouped by chromosome and window
# The table structure contains:
#   - chr: chromosome number
#   - window: 100kb genomic window identifier
#   - Precombinations: paternal recombination events
#   - Mrecombinations: maternal recombination events
#   - recombinations: total recombination events
#   - *_female/*_male: sex-specific recombination counts
#--------------------------------------------------------------------

# SQL query to create f2inheritance_100k_recombinations table
sql_create_recombination_table = "
CREATE TABLE f2inheritance_100k_recombinations AS
SELECT chr, window,
  SUM(CASE WHEN Pallelecount > 1 THEN 1 ELSE 0 END) AS Precombinations,
  SUM(CASE WHEN Mallelecount > 1 THEN 1 ELSE 0 END) AS Mrecombinations,
  SUM(CASE WHEN Mallelecount > 1 THEN 1 ELSE 0 END) + 
  SUM(CASE WHEN Pallelecount > 1 THEN 1 ELSE 0 END) AS recombinations,
  SUM(CASE WHEN Pallelecount > 1 AND f2 % 2 = 0 THEN 1 ELSE 0 END) AS Precombinations_female,
  SUM(CASE WHEN Mallelecount > 1 AND f2 % 2 = 0 THEN 1 ELSE 0 END) AS Mrecombinations_female,
  SUM(CASE WHEN Mallelecount > 1 AND f2 % 2 = 0 THEN 1 ELSE 0 END) + 
  SUM(CASE WHEN Pallelecount > 1 AND f2 % 2 = 0 THEN 1 ELSE 0 END) AS recombinations_female,
  SUM(CASE WHEN Pallelecount > 1 AND f2 % 2 = 1 THEN 1 ELSE 0 END) AS Precombinations_male,
  SUM(CASE WHEN Mallelecount > 1 AND f2 % 2 = 1 THEN 1 ELSE 0 END) AS Mrecombinations_male,
  SUM(CASE WHEN Mallelecount > 1 AND f2 % 2 = 1 THEN 1 ELSE 0 END) + 
  SUM(CASE WHEN Pallelecount > 1 AND f2 % 2 = 1 THEN 1 ELSE 0 END) AS recombinations_male
FROM f2inheritance_100k_allelecount
GROUP BY chr, window"

# Execute the query to create recombination statistics table
cat("Creating f2inheritance_100k_recombinations table...\n")
result <- dbSendQuery(con, sql_create_recombination_table)
dbClearResult(result)

#--------------------------------------------------------------------
# Create window-based recombination statistics
# This analysis identifies recombination events based on fragment inheritance patterns
# Step 1: Create temporary table to mark recombination status
# Step 2: Generate final recombination statistics
# Step 3: Add performance indexes
#--------------------------------------------------------------------

# Step 1: Create temporary table to mark recombination status
# Individuals with more than 2 fragments in a window are considered recombined
sql_temp_recombination = "
CREATE TEMPORARY TABLE temp_recombination_status AS
SELECT 
    f2, 
    chr, 
    window,
    CASE WHEN COUNT(*) > 2 THEN 1 ELSE 0 END AS is_recombined
FROM 
    f2fragmentinheritance_window100k
GROUP BY 
    f2, chr, window"

cat("Creating temporary recombination status table...\n")
result <- dbSendQuery(con, sql_temp_recombination)
dbClearResult(result)

# Step 2: Create final recombination statistics table
# This table contains aggregated statistics for each genomic window
sql_final_recombination = "
CREATE TABLE window_recombination_statistics AS
SELECT 
    chr, 
    window,
    SUM(is_recombined) AS recombined_count,
    COUNT(*) - SUM(is_recombined) AS non_recombined_count,
    COUNT(*) AS total_individuals,
    ROUND(SUM(is_recombined) / COUNT(*) * 100, 2) AS recombination_rate
FROM 
    temp_recombination_status
GROUP BY 
    chr, window
ORDER BY 
    chr, window"

cat("Creating window recombination statistics table...\n")
result <- dbSendQuery(con, sql_final_recombination)
dbClearResult(result)

# Step 3: Add index to optimize queries
# Index on chromosome and window for faster data retrieval
sql_create_index = "CREATE INDEX idx_window_recomb_stats ON window_recombination_statistics(chr, window)"

cat("Creating performance index...\n")
result <- dbSendQuery(con, sql_create_index)
dbClearResult(result)

# Step 4: Drop temporary table to clean up
sql_drop_temp = "DROP TABLE temp_recombination_status"

cat("Cleaning up temporary tables...\n")
result <- dbSendQuery(con, sql_drop_temp)
dbClearResult(result)

#--------------------------------------------------------------------
# Verification and summary statistics
#--------------------------------------------------------------------

# Check if tables were created successfully
cat("Verifying table creation...\n")
tables_query = "SHOW TABLES LIKE '%recombination%'"
tables_result = dbGetQuery(con, tables_query)
print(tables_result)

# Get summary statistics from the created tables
cat("\nSummary of f2inheritance_100k_recombinations:\n")
summary_query1 = "SELECT COUNT(*) as total_windows, 
                         AVG(recombinations) as avg_recombinations,
                         MAX(recombinations) as max_recombinations,
                         MIN(recombinations) as min_recombinations
                  FROM f2inheritance_100k_recombinations"
summary1 = dbGetQuery(con, summary_query1)
print(summary1)

cat("\nSummary of window_recombination_statistics:\n")
summary_query2 = "SELECT COUNT(*) as total_windows,
                         AVG(recombination_rate) as avg_recombination_rate,
                         MAX(recombination_rate) as max_recombination_rate,
                         MIN(recombination_rate) as min_recombination_rate
                  FROM window_recombination_statistics"
summary2 = dbGetQuery(con, summary_query2)
print(summary2)

# Close database connection
dbDisconnect(con)