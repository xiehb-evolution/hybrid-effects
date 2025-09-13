require(RMySQL)
require(ggplot2)
require(sqldf)
require(ggpubr)

sql="SELECT concat(\"chr\",chr) as chr,(pos1+pos2)/100000 as winid,class,count(*) as count
FROM f2recombination where chr<23 and class=\"P\"
group by chr,(pos1+pos2)/100000,class"

result = dbGetQuery(con, sql)

ggplot(data = result, aes(x=winid,y=count)) + geom_point(cex=1) + facet_wrap(~chr, ncol = 6,scales = "free_x")


sql="SELECT chr as chrnew,(pos1+pos2)/100000 as winid,class,count(*) as count
FROM f2recombination where chr<23 and class=\"P\"
group by chr,(pos1+pos2)/100000,class
order by chr"

result = dbGetQuery(con, sql)
plotA = ggplot(data = result, aes(x=winid,y=count)) + 
  geom_point(cex=1) + facet_wrap(~chrnew, ncol = 6,scales = "free_x") + 
  ylab("Number of recombination events in paternal genome") + 
  xlab("Chromosomal position (Mb)") + 
  scale_x_continuous(labels = ~ .x / 10)

sql="SELECT chr as chrnew,(pos1+pos2)/100000 as winid,class,count(*) as count
FROM f2recombination where chr<23 and class=\"M\"
group by chr,(pos1+pos2)/100000,class"
result = dbGetQuery(con, sql)
plotB = ggplot(data = result, aes(x=winid,y=count)) + 
  geom_point(cex=1) + facet_wrap(~chrnew, ncol = 6,scales = "free_x") +
  ylab("Number of recombination events in maternal genome") + 
  xlab("Chromosomal position (Mb)") + 
  scale_x_continuous(labels = ~ .x / 10)

ggarrange(plotA,plotB,ncol=1,labels=c("A","B"))





