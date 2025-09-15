setwd("E:\\xiehb_sync\\F2\\data\\sscrofa11.data")

require(sqldf)

#load pedigree data
haplotype = read.table("chr1.phased.duohmm.haps")
sample = read.table("chr1.phased.duohmm.sample",skip=2)
f2id = as.data.frame(as.matrix(sample[77:654,2]))

# The following data structure saves the family information
# the column values are the sample id
#-------------------------------------------------------------------------  
# table format (7 columns):
#-------------------------------------------------------------------------  
#       F2   F1male   F0male   F0female   F1female  F0male  F0female
#     ^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^
#     f2id          family1                       family2
#-------------------------------------------------------------------------  

# -------- It is much safer to use the following code to generate pedigree information
f1father = sqldf("select a.V1,b.V4 from f2id a, sample b where a.V1=b.V2",drv="SQLite")
f1mother = sqldf("select a.V1,b.V5 from f2id a, sample b where a.V1=b.V2",drv="SQLite")
family1 = sqldf("select a.V1,b.V2,b.V4,b.V5 from f1father a, sample b where a.V4=b.V2",drv="SQLite")
family2 = sqldf("select a.V1,b.V2,b.V4,b.V5 from f1mother a, sample b where a.V5=b.V2",drv="SQLite")
pedigree = sqldf("select a.V1,a.V2,a.V4,a.V5,b.V2,b.V4,b.V5 from family1 a,family2 b where a.V1=b.V1",drv="SQLite")

pedigree_tmp = as.data.frame(matrix(as.vector(as.matrix(pedigree)),byrow=T,ncol=1))
allsample = as.vector(sample[[2]])
allsample = as.data.frame(cbind(1:length(allsample),allsample))
allsamplepos = sqldf("select a.V1 from pedigree_tmp b left join allsample a on b.V1 = a.allsample",drv="SQLite")
allsamplepos = as.vector(as.matrix(allsamplepos))

#pedigree_info saves the column # for each individual
pedigree_info = matrix(allsamplepos,byrow=F,nrow=nrow(pedigree))

# pedigree_haps saves the column # for chromosomes in the haplotype (data frame)
pedigree_haps1 = (pedigree_info-1)*2+6 
pedigree_haps2 = (pedigree_info-1)*2+7
pedigree_haps = cbind(pedigree_haps1,pedigree_haps2)[,c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]

# add SNP coordinates column # to pedigree_haps 
pedigree_haps = cbind(3,pedigree_haps)
pedigree_haps = cbind(1,pedigree_haps)

#determine the difference between two haplotypes
# m and n are numbers indicating the column # in the data structure
GetScore = function(m,n)
{
  return(abs(pedigree_haplotype[[m]]-pedigree_haplotype[[n]]))
}

# To figure out which part of the haplotype in a chromosome (x) is inherited from which part of the haplotypes of parental chromosomes y and z.
# The primary role of this function is to identify putative recombination events
# The return value is the haplotype of x encoded by y and z
# x, y, z is a column # specifying different haplotypes in the data frame for a specific family
StepHap = function(x,y,z,n)
{
  mismatchY = GetScore(x,y)
  mismatchZ = GetScore(x,z)
  
  finalseq = 0
  localsites = n
  
  finalseq=rep(0,localsites)
  lastparent = 0
  lastblockend = 0
  
  while(length(finalseq[finalseq==0])>0)
  {
    vectorY = which(mismatchY[(lastblockend+1):length(mismatchY)]==1)
    vectorZ = which(mismatchZ[(lastblockend+1):length(mismatchZ)]==1)
    
    firstY = min(vectorY) + lastblockend
    firstZ = min(vectorZ) + lastblockend
    firstY
    firstZ
    
    if(is.infinite(firstY) && is.infinite(firstZ)==F)
    {
      finalseq[(lastblockend+1):length(mismatchY)] = y
      break
    }
    if(is.infinite(firstZ) && is.infinite(firstY)==F)
    {
      finalseq[(lastblockend+1):length(mismatchZ)] = z
      break
    }
    if(is.infinite(firstY) && is.infinite(firstZ))
    {
      finalseq[(lastblockend+1):length(mismatchZ)] = lastparent
      break
    }
    
    if(firstY==firstZ)
    {
      if(lastparent>0)
      {
        finalseq[(lastblockend+1):firstY] = lastparent
      }
      else
      {
        finalseq[(lastblockend+1):firstY] = y
      }
      lastblockend = firstY
    } else
    {
      if(firstZ > firstY)
      {
        finalseq[(lastblockend+1):(firstZ-1)] = z
        lastblockend = firstZ - 1
        lastparent = z
      } else
      {
        finalseq[(lastblockend+1):(firstY-1)] = y
        lastblockend = firstY - 1
        lastparent = y
      }
    }
  }
  
  return(finalseq)
}

#--------------------------------------------------------------------
# This function is to correct genotyping errors & switching errors
# It takes three steps
correcting_haplotype = function(onehap,columnflag,minfraglen,mindiffs)
{
  alleles = unique(onehap)
  
  #step1
  #remove genotyping errors
  if(length(onehap[onehap==alleles[1]])<=3)
  {
    onehap[onehap==alleles[1]] = alleles[2]
    return(onehap)
  }
  if(length(onehap[onehap==alleles[2]])<=3)
  {
    onehap[onehap==alleles[2]] = alleles[1]
    return(onehap)
  }
  
  #step2
  #remove switching errors with short genomic spans
  finderr=1
  i=1:(length(onehap)-1)
  j=i+1
  while(finderr==1)
  {
    finderr=0
    seps = c(0,which(onehap[i]!=onehap[j]),length(onehap))
    seps
    if(length(seps)>=3)
    {
      frags = cbind(seps[1:(length(seps)-1)]+1,seps[2:length(seps)])
      fraglen = frags[,2]-frags[,1]+1
      shortspanid = which(fraglen<=minfraglen)
      if(length(shortspanid)>0)
      {
        for(k in 1:length(shortspanid))
        {
          finderr=1
          onehap[frags[shortspanid[k],1]:frags[shortspanid[k],2]] = columnflag - onehap[frags[shortspanid[k],1]:frags[shortspanid[k],2]]
        }
      }
    }
  }
  
  #step3
  #remove switching-cross with a few informative markers
  seps = c(0,which(onehap[i]!=onehap[j]),length(onehap))
  if(length(seps)>=3)
  {
    frags = cbind(seps[1:(length(seps)-1)]+1,seps[2:length(seps)])
    fraglen = frags[,2]-frags[,1]+1
    hapdiff = pedigree_haplotype[[(columnflag-1)/2]]-pedigree_haplotype[[(columnflag+1)/2]]
    for(k in 1:nrow(frags))
    {
      fragdiff = sum(abs(hapdiff[frags[k,1]:frags[k,2]]))
      if(fragdiff<=mindiffs)
      {
        onehap[frags[k,1]:frags[k,2]] = columnflag - onehap[frags[k,1]:frags[k,2]]
      }
    }
  }
  return(onehap)
}

#for debug
#check family for a specific F2 individual, row # used here
# pedigreeid=0
# pedigreeid=pedigreeid+1
# pedigree_haplotype = haplotype[,pedigree_haps[pedigreeid,]]
# score1 = sum(GetScore(3,5)+GetScore(4,11))
# score2 = sum(GetScore(3,5)+GetScore(4,12))
# score3 = sum(GetScore(3,6)+GetScore(4,11))
# score4 = sum(GetScore(3,6)+GetScore(4,12))
# score5 = sum(GetScore(4,5)+GetScore(3,11))
# score6 = sum(GetScore(4,5)+GetScore(3,12))
# score7 = sum(GetScore(4,6)+GetScore(3,11))
# score8 = sum(GetScore(4,6)+GetScore(3,12))
# score = c(score1,score2,score3,score4,score5,score6,score7,score8)
# optpos = which(score==min(score))
# optpos
# 
# chrhap1=0
# chrhap2=0
# if(optpos<=4)
# {
#   x=3
#   y=5
#   z=6
#   chrhap1 = StepHap(x,y,z,nrow(haplotype))
#   x=4
#   y=11
#   z=12
#   chrhap2 = StepHap(x,y,z,nrow(haplotype))
# } else
# {
#   x=4
#   y=5
#   z=6
#   chrhap1 = StepHap(x,y,z,nrow(haplotype))
#   x=3
#   y=11
#   z=12
#   chrhap2 = StepHap(x,y,z,nrow(haplotype))
# }
# 
# plot(correcting_haplotype(chrhap1,11,100,10),main=pedigreeid)
# plot(chrhap1,main=pedigreeid)
# plot(chrhap2,main=pedigreeid)
# pedigree_info[pedigreeid,]
# pedigree[pedigreeid,]


#--------------------------------------------------------
#   Constructing Haplotype
#--------------------------------------------------------
# Recode F2 haplotype using F1 haplotypes, j is the chromosome ID (for the pig, it ranges from 1 to 18)
j = 1
for(j in 1:18)
{
    haplotype = read.table(sprintf("chr%d.phased.duohmm.haps",j))
    sample = read.table(sprintf("chr%d.phased.duohmm.sample",j),skip=2)

    F2RecodedByF1 = haplotype
    for(pedigreeid in 1:nrow(f2id))
    {
      cat(sprintf("pedigreeid=%d\n",pedigreeid))
      pedigree_haplotype = haplotype[,pedigree_haps[pedigreeid,]]
      score1 = sum(GetScore(3,5)+GetScore(4,11))
      score2 = sum(GetScore(3,5)+GetScore(4,12))
      score3 = sum(GetScore(3,6)+GetScore(4,11))
      score4 = sum(GetScore(3,6)+GetScore(4,12))
      score5 = sum(GetScore(4,5)+GetScore(3,11))
      score6 = sum(GetScore(4,5)+GetScore(3,12))
      score7 = sum(GetScore(4,6)+GetScore(3,11))
      score8 = sum(GetScore(4,6)+GetScore(3,12))
      
      score = c(score1,score2,score3,score4,score5,score6,score7,score8)
      score
      optpos = which(score==min(score))
      optpos = optpos[1]
      optpos
      chrhap1=0
      chrhap2=0
      
      if(optpos<=4)
      {
        x=3
        y=5
        z=6
        chrhap1 = StepHap(x,y,z,nrow(haplotype))
        x=4
        y=11
        z=12
        chrhap2 = StepHap(x,y,z,nrow(haplotype))
      } else
      {
        x=4
        y=5
        z=6
        chrhap1 = StepHap(x,y,z,nrow(haplotype))
        x=3
        y=11
        z=12
        chrhap2 = StepHap(x,y,z,nrow(haplotype))
      }
      chrhap2
      F2RecodedByF1[,pedigree_haps[pedigreeid,][3]] = correcting_haplotype(chrhap1,11,50,10)
      F2RecodedByF1[,pedigree_haps[pedigreeid,][4]] = correcting_haplotype(chrhap2,23,50,10)
    }

    #Save the data to f2inheritance.txt
    for(pedigreeid in 1:nrow(f2id))
    {
      x = cbind(f2id[pedigreeid,1],F2RecodedByF1[,c(1,3)],F2RecodedByF1[,pedigree_haps[pedigreeid,][3]], F2RecodedByF1[,pedigree_haps[pedigreeid,][4]])
      write.table(x, "f2inheritance.txt",append=T,quote=F,row.names=F,col.names=F,sep="\t")
    }

    # Determine the recombination breakpoints
    for(pedigreeid in 1:nrow(f2id))
    {
      faminfo = pedigree[pedigreeid,c(1,2,5)]
      
      paternal = F2RecodedByF1[,pedigree_haps[pedigreeid,][3]]
      maternal = F2RecodedByF1[,pedigree_haps[pedigreeid,][4]]
      paternal1 = paternal[-length(paternal)]
      paternal2 = paternal[-1]
      diffpos = which(paternal1!=paternal2)
      if(length(diffpos)>0)
      {
        for (k in 1:length(diffpos))
        {
          x = cbind(faminfo,F2RecodedByF1[diffpos[k],1],F2RecodedByF1[diffpos[k],3], F2RecodedByF1[diffpos[k]+1,3],"P")
          write.table(x, "f2recombination.txt",append=T,quote=F,row.names=F,col.names=F,sep="\t")
        }
      }
      maternal1 = maternal[-length(maternal)]
      maternal2 = maternal[-1]
      diffpos = which(maternal1!=maternal2)
      if(length(diffpos)>0)
      {
        for (k in 1:length(diffpos))
        {
          x = cbind(faminfo,F2RecodedByF1[diffpos[k],1],F2RecodedByF1[diffpos[k],3], F2RecodedByF1[diffpos[k]+1,3],"M")
          write.table(x, "f2recombination.txt",append=T,quote=F,row.names=F,col.names=F,sep="\t")
        }
      }
    }

}

#Save the f2inheritance.txt to a MySQL server for further analysis
f2inheritance = read.table("f2inheritance.txt", head=F)
names(f2inheritance) = c("f2", "chr", "pos", "f1father", "f1mother")
require(RMySQL)
mysqlserver="10.0.16.70"
mysqlport=3306
mysqldbname="speciation"
mysqluser="xiehb"
mysqlpassword="password"
con<-dbConnect(dbDriver("MySQL"),host=mysqlserver,port=mysqlport,dbname=mysqldbname,user=mysqluser,password=mysqlpassword)

dbWriteTable(con, "f2inheritance", f2inheritance,overwrite = TRUE, row.names = FALSE)
sql = "create index idx1 on f2inheritance(f2, f1father, f1mother)"
ret = dbSendQuery(con, sql)
sql = "create index idx2 on f2inheritance(chr, pos)"
ret = dbSendQuery(con, sql)
