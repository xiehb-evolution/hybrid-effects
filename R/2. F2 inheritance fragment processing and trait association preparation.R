require(RMySQL)
require(sqldf)
require(ggplot2)

#--------------------------------------------------------
# Generate the f2fragmentinheritance from f2inheritance
# This step is to identify genomic fragments that are inherited as a whole
# The paternal and maternal genomes are calculated independently
# The resulting file is f2inheritance.length.txt and a relevant table is f2fragmentinheritance in database
#--------------------------------------------------------

# Get all unique F2 individual IDs from chromosome 1 data
sql=sprintf("select distinct f2 from f2inheritance where chr=1")
result = dbGetQuery(con, sql)
f2all  = result[[1]]

# Process each chromosome (1-18) for each F2 individual
for(chr in c(1:18))
{
  for(f2id in 1:length(f2all))
  {
      f2=f2all[f2id]
      
      # Extract inheritance data for current F2 individual and chromosome, ordered by position
      sql=sprintf("select * from f2inheritance where chr=%d and f2=%d order by pos",chr,f2)
      result = dbGetQuery(con, sql)
      
      # Identify maternal inheritance blocks
      # Find positions where maternal inheritance changes between adjacent SNPs
      d1 = result$f1mother[-length(result$f1mother)]  # all positions except last
      d2 = result$f1mother[-1]                        # all positions except first
      end = which(d1!=d2)                             # positions where inheritance changes
      start = c(1,end+1)                              # start positions of inheritance blocks
      end = c(end,length(result$f1mother))            # end positions of inheritance blocks
      # Create maternal inheritance fragments: F2_ID, chromosome, start_pos, end_pos, inheritance_pattern, origin
      m = cbind(f2,chr,result[start,3],result[end,3],result[start,5]%%2,'M')
      
      # Identify paternal inheritance blocks using same approach as maternal
      d1 = result$f1father[-length(result$f1father)]
      d2 = result$f1father[-1]
      end = which(d1!=d2)
      start = c(1,end+1)
      end = c(end,length(result$f1father))
      # Create paternal inheritance fragments
      p = cbind(f2,chr,result[start,3],result[end,3],result[start,4]%%2,'P')
      
      # Combine maternal and paternal inheritance data
      data = rbind(m,p)
      write.table(data,"f2inheritance.length.txt",sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
}

# Load the inheritance fragment data and create database table
outval = read.table("f2inheritance.length.txt",head=F)
names(outval)=c("f2","chr","start","end","inheritance","origin")
dbRemoveTable(con,"f2fragmentinheritance")
dbWriteTable(con,"f2fragmentinheritance",outval,append=T,row.names=F)

# Create database indexes for efficient querying
sql = "create index idx1 on f2fragmentinheritance (inheritance, origin);"
ret = dbSendQuery(con, sql)
sql = "create index idx2 on f2fragmentinheritance (f2, chr, start, end, inheritance, origin);"
ret = dbSendQuery(con, sql)
sql = "create index idx3 on f2fragmentinheritance (chr);"
ret = dbSendQuery(con, sql)



#--------------------------------------------------------
# Generating f2fragmentinheritance_window100k table in mysql server
# from the table f2fragmentinheritance
# This step generates the genotypes of 100-kb sliding windows for each individual
# Convert genomic positions to 100kb window coordinates for analysis
#--------------------------------------------------------

# Convert start and end positions to 100kb window numbers
sql = "select f2,chr,floor(start/100000) as startwin,floor(end/100000) as endwin,inheritance,origin from f2fragmentinheritance;"
result = dbGetQuery(con, sql)
result[1:10,]

# Set working directory for output files
#setwd("E:\\xiehb_sync\\F2\\data\\sscrofa11.data")
setwd("C:\\Users\\xiehb\\xiehb_sync\\F2\\data\\sscrofa11.data")

# For each inheritance fragment, create entries for all 100kb windows it spans
for(i in 1:nrow(result))
{
  # Generate window entries from startwin to endwin for each fragment
  write.table(cbind(result[i,1:2],result[i,3]:result[i,4],result[i,5:6]),"f2fragmentinheritance_window100k.txt",quote=F,row.names=F,col.names=F,append=T)
}

# Load window data and create database table
outval = c()
outval = read.table("f2fragmentinheritance_window100k.txt",head=F)
names(outval)=c("f2","chr","window","inheritance","origin")
dbRemoveTable(con,"f2fragmentinheritance_window100k")
dbWriteTable(con,"f2fragmentinheritance_window100k",outval,append=T,row.names=F)

# Optimize column data types for storage efficiency
sql = "ALTER TABLE f2fragmentinheritance_window100k
MODIFY COLUMN f2 int NULL DEFAULT NULL FIRST,
MODIFY COLUMN chr tinyint NULL DEFAULT NULL AFTER f2,
MODIFY COLUMN window smallint NULL DEFAULT NULL AFTER chr,
MODIFY COLUMN inheritance tinyint NULL DEFAULT NULL AFTER window,
MODIFY COLUMN origin char(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci NULL DEFAULT NULL AFTER inheritance"
ret = dbSendQuery(con, sql)

# Create index for efficient querying
sql = "create index idx1 on f2fragmentinheritance_window100k (f2, chr, window, inheritance);"
ret = dbSendQuery(con, sql)

# Remove sex chromosome data for males (f2%2=1 indicates males, chr=23 is sex chromosome)
sql = "delete from f2fragmentinheritance_window100k where f2%2=1 and chr=23 and origin=\"P\""
ret = dbSendQuery(con, sql)


#--------------------------------------------------------
# Fill missing values in f2fragmentinheritance_window100k
# There are missing values due to some exceptions where there are no informative SNPs in some sliding windows
# The missing values are filled with genotypes of flanking windows
#--------------------------------------------------------

# Get all F2 individuals
sql = "select distinct f2 from f2fragmentinheritance_window100k"
f2 = dbGetQuery(con, sql)

# Process each F2 individual for each autosomal chromosome
for(i in 1:nrow(f2))
{
  for(j in 1:18)  # Only autosomal chromosomes
  {
    # Get paternal inheritance windows for current F2 and chromosome
    sql = sprintf("select window,inheritance from f2fragmentinheritance_window100k 
                  where f2=%s and chr=%d and origin='P' group by window,inheritance order by window",f2[i,1],j)
    windowdata = dbGetQuery(con, sql)
    window = as.vector(as.matrix(windowdata[,1]))
    inheritance = as.vector(as.matrix(windowdata[,2]))
    
    # Find missing windows (gaps larger than 1 in window sequence)
    misspos = which( (window[-1] - window[-length(window)])!=1 & (window[-1] - window[-length(window)])!=0)
    misspos
    
    # Fill missing windows for paternal inheritance
    if(length(misspos)>0)
    {
      print(sprintf("f2=%d chr=%d misspos=%d P",f2[i,1],j,misspos[1]))
      for(x in 1:length(misspos))
      {
        k = misspos[x]
        # Calculate recombination position as midpoint between flanking windows
        recpos = floor((window[k+1]+window[k])/2)
        
        # Fill windows before recombination position with previous inheritance pattern
        if((window[k]+1)<recpos)
        {
          for(m in (window[k]+1):(recpos-1))
          {
            sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'P')",f2[i,1],j,m,inheritance[k])
            dbSendQuery(con, sql)
            print(sprintf("%s chr=%d window=%d 'P'",f2[i,1],j,m))
          }
        }
        
        # Insert recombination position with both inheritance patterns
        sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'P')",f2[i,1],j,recpos,inheritance[k])
        dbSendQuery(con, sql)
        print(sprintf("%s chr=%d window=%d 'P' recpos1",f2[i,1],j,recpos))
        
        sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'P')",f2[i,1],j,recpos,1-inheritance[k])
        dbSendQuery(con, sql)
        print(sprintf("%s chr=%d window=%d 'P' recpos2",f2[i,1],j,recpos))
        
        # Fill windows after recombination position with next inheritance pattern
        if((window[k+1]-1)>recpos)
        {
          for(m in (recpos+1):(window[k+1]-1))
          {
            sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'P')",f2[i,1],j,m,inheritance[k+1])
            dbSendQuery(con, sql)
            print(sprintf("%s chr=%d window=%d 'P'",f2[i,1],j,m))
          }
        }        
      }
    }

    # Repeat the same process for maternal inheritance
    sql = sprintf("select window,inheritance from f2fragmentinheritance_window100k 
                  where f2=%s and chr=%d and origin='M' group by window,inheritance order by window",f2[i,1],j)
    windowdata = dbGetQuery(con, sql)
    window = as.vector(as.matrix(windowdata[,1]))
    inheritance = as.vector(as.matrix(windowdata[,2]))
    misspos = which( (window[-1] - window[-length(window)])!=1 & (window[-1] - window[-length(window)])!=0)
    misspos
    
    # Fill missing windows for maternal inheritance (same logic as paternal)
    if(length(misspos)>0)
    {
      print(sprintf("f2=%d chr=%d misspos=%d M",f2[i,1],j,misspos[1]))
      for(x in 1:length(misspos))
      {
        k = misspos[x]
        recpos = floor((window[k+1]+window[k])/2)
        if((window[k]+1)<recpos)
        {
          for(m in (window[k]+1):(recpos-1))
          {
            sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'M')",f2[i,1],j,m,inheritance[k])
            dbSendQuery(con, sql)
            print(sprintf("%s chr=%d window=%d 'M'",f2[i,1],j,m))
          }
        }
        sql = sprintf("insert into f2fragmentinheritance_window100k 
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'M')",f2[i,1],j,recpos,inheritance[k])
        dbSendQuery(con, sql)
        print(sprintf("%s chr=%d window=%d 'M' recpos1",f2[i,1],j,recpos))
        
        sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'M')",f2[i,1],j,recpos,1-inheritance[k])
        dbSendQuery(con, sql)
        print(sprintf("%s chr=%d window=%d 'M' recpos2",f2[i,1],j,recpos))
        
        if((window[k+1]-1)>recpos)
        {
          for(m in (recpos+1):(window[k+1]-1))
          {
            sql = sprintf("insert into f2fragmentinheritance_window100k
                    (f2,chr,window,inheritance,origin)
                    values (%s,%d,%d,%d,'M')",f2[i,1],j,m,inheritance[k+1])
            dbSendQuery(con, sql)
            print(sprintf("%s chr=%d window=%d 'M'",f2[i,1],j,m))
          }
        }        
      }
    }   
        
  }
}



#------------------------------------------------------------------------------------------------------------
# Remove windows with recombination breakpoints
# The table f2fragmentinheritance_window100k contains sliding windows on which some F2 carry recombination breakpoints
# The recombination breaks within a window will cause the genotype of the window not unique
# The samples should be removed on such windows in the f2fragmentinheritance_window table
# The resulting table is f2fragmentinheritance_window100k_without_recombination
#------------------------------------------------------------------------------------------------------------

# Create a copy of the window table
sql = "create table f2fragmentinheritance_window100k_without_recombination select * from f2fragmentinheritance_window100k"
dbSendQuery(con, sql)

# Create index for efficient deletion operations
sql = "create index idx1 on f2fragmentinheritance_window100k_without_recombination(f2, chr, window, inheritance)"
dbSendQuery(con, sql)

# Identify windows with recombination breakpoints (more than 2 entries per F2-chr-window combination)
# Normal windows should have exactly 2 entries (one paternal, one maternal)
sql = "select f2,chr,window from f2fragmentinheritance_window100k group by f2,chr,window having count(*)>2"
result = dbGetQuery(con, sql)

# Remove all entries for windows containing recombination breakpoints
for(i in 1:nrow(result))
{
  sql = sprintf("delete from f2fragmentinheritance_window100k_without_recombination 
                  where f2=%d and chr=%d and window=%d",result[i,1],result[i,2],result[i,3])
  dbSendQuery(con, sql)
}


#------------------------------------------------------------------------------------------------------------
# Convert inheritance data to combined paternal-maternal format
# Convert f2fragmentinheritance_window100k_without_recombination into 
# f2fragmentinheritance_window100k_without_recombination_PM for subsequent analysis
# This creates a single row per F2-chr-window with both paternal and maternal inheritance patterns
#------------------------------------------------------------------------------------------------------------
sql = "create table f2fragmentinheritance_window100k_without_recombination_PM
SELECT a.f2,a.chr,a.window,a.inheritance as paternalinheritance,b.inheritance as maternalinheritance
FROM f2fragmentinheritance_window100k_without_recombination a,f2fragmentinheritance_window100k_without_recombination b
where a.chr=b.chr and a.f2=b.f2 and a.window=b.window and a.origin=\"P\" and b.origin=\"M\"
order by f2,chr,window"
dbSendQuery(con, sql)



#------------------------------------------------------------------------------------------------------------
# Combine inheritance and trait data
# Create a table that combines inheritance pattern data with trait information for each F2 individual
# This joins fragment inheritance data with trait measurements using the F2 individual ID as the key
#------------------------------------------------------------------------------------------------------------
sql = "create table window100k_single_site_trait
SELECT a.*,b.trait_id,b.trait_name,b.trait_value
FROM f2fragmentinheritance_window100k_without_recombination_PM a, f2_trait_name_trait_value_copy1 b
where a.f2=b.f2"
dbSendQuery(con, sql)

#------------------------------------------------------------------------------------------------------------
# Calculate trait statistics by inheritance pattern
# Create a table calculating sex-specific trait statistics (counts and averages) grouped by trait, chromosome, 
# window, and inheritance patterns (paternal/maternal)
# This enables QTL analysis by comparing trait values across different inheritance patterns
#------------------------------------------------------------------------------------------------------------
sql = "create table window100k_single_site_trait_stat
SELECT trait_id,chr,window,paternalinheritance,maternalinheritance,
sum(case f2%2=1 when 1 then 1 else 0 end) as malecount,
sum(case f2%2=0 when 1 then 1 else 0 end) as femalecount,
sum(case f2%2=1 when 1 then 1 else 0 end * trait_value)/sum(case f2%2=1 when 1 then 1 else 0 end) as maletrait,
sum(case f2%2=0 when 1 then 1 else 0 end * trait_value)/sum(case f2%2=0 when 1 then 1 else 0 end) as femaletrait
FROM window100k_single_site_trait
group by trait_id,chr,window,paternalinheritance,maternalinheritance"
dbSendQuery(con, sql)

# Create indexes to optimize queries on the statistical table
sql1 = "CREATE INDEX idx1 ON window100k_single_site_trait_stat(trait_id, chr, window, paternalinheritance, maternalinheritance);"
sql2 = "CREATE INDEX idx2 ON window100k_single_site_trait_stat(trait_id, chr, window, maternalinheritance, paternalinheritance);"
dbSendQuery(con, sql1)
dbSendQuery(con, sql2)


#------------------------------------------------------------------------------------------------------------
# Analysis of Trait Variation by Inheritance Pattern
# This section calculates the deviation in mean trait values between different inheritance patterns.
# The analysis includes sex-specific differences and identifies genomic regions
# where changes in inheritance patterns are associated with trait variation.
#------------------------------------------------------------------------------------------------------------
sql = "create table window100k_single_site_trait_stat_mutant_deviation_from_mean_all
SELECT a.*,b.paternalinheritance as paternalinheritance_mutant,b.maternalinheritance as maternalinheritance_mutant,
b.malecount-a.malecount as malecount_diff,
b.femalecount-a.femalecount as femalecount_diff,
b.maletrait-a.maletrait as maletrait_diff,
b.femaletrait-a.femaletrait as femaletrait_diff,
(a.maletrait-c.traitmean)/c.traitmean as maledev,
(a.femaletrait-d.traitmean)/d.traitmean as femaledev,
(b.maletrait-c.traitmean)/c.traitmean as malemutantdev,
(b.femaletrait-d.traitmean)/d.traitmean as femalemutantdev
FROM window100k_single_site_trait_stat a,window100k_single_site_trait_stat b, f2_trait_name_trait_value_copy1_mean c, f2_trait_name_trait_value_copy1_mean d
where a.trait_id=b.trait_id and a.chr=b.chr and a.window=b.window and a.trait_id=c.trait_id and c.sex=1 and a.trait_id=d.trait_id and d.sex=0 and not (a.paternalinheritance=b.paternalinheritance and a.maternalinheritance=b.maternalinheritance)"
dbSendQuery(con, sql)


# Filtering for Single-Locus Effects
# This step filters the results to include only autosomal chromosomes and inheritance pattern
# pairs that differ by exactly one allele. This focuses the analysis on the effects of
# single-locus changes rather than on more complex interactions.
sql = "create table window100k_single_site_trait_stat_mutant_deviation_from_mean
select * from window100k_single_site_trait_stat_mutant_deviation_from_mean_all
where abs(paternalinheritance+maternalinheritance-paternalinheritance_mutant-maternalinheritance_mutant)<2 and 
abs(paternalinheritance+maternalinheritance-paternalinheritance_mutant-maternalinheritance_mutant)>0 and chr<23;"
dbSendQuery(con, sql)

# Create indexes to optimize queries on the deviation analysis table
sql1 = "create index idx1 on window100k_single_site_trait_stat_mutant_deviation_from_mean(trait_id,chr,window,paternalinheritance,maternalinheritance,paternalinheritance_mutant, maternalinheritance_mutant);"
sql2 = "create index idx2 on window100k_single_site_trait_stat_mutant_deviation_from_mean(chr,window,paternalinheritance,maternalinheritance,paternalinheritance_mutant,maternalinheritance_mutant,trait_id);"
dbSendQuery(con, sql1)
dbSendQuery(con, sql2)

