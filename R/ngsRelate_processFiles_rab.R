### do ngsRelate values (r>0.8) show up as clones? 
library(fBasics)
rm(list=ls())
### meta of adults####
inds <- readLines("data/268adults99seeds_clean")
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2),
                  grid=substr(inds,3,3),
                  ind2=substr(inds,5,6))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "cur"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "nil"
env$site.nice[env$site.nice=="W"] <- "wes"
env$site.nice[env$site.nice=="L"] <- "lyn"
env$site.nice <- factor(env$site.nice)
env$quad.name <- paste(env$site.nice,".",env$SvD,env$grid,"_",env$ind2,sep="")
env$site.nice.depth <- paste(env$site.nice,env$SvD,sep="_")
env$site.nice.depth <- factor(env$site.nice.depth)
env$site.nice.depth <- factor(env$site.nice.depth,levels(env$site.nice.depth)[c(1:4,7,8,5,6)])

### ngsRelate data
filename = paste("data/ngsRelateFiles/",c(
  "268adults99seeds95perc_newres.gz",
  "268adults99seeds90perc_newres.gz",
  "268adults99seeds80perc_newres.gz",
  "268adults99seeds50perc_newres.gz"
  ),sep="")

out = c()
pdf("output/histogram_ngsRelate_rab.pdf")
for(i in 1:length(filename))
{
rab <- read.delim(file=filename[i])[,c("a","b","rab")]
rab$envA = env[rab$a+1,]
rab$envB = env[rab$b+1,]

### relatedness between adults and seeds and between seeds and seeds 
rab$adult.seed2 = paste(rab$envA$adult.seed,rab$envB$adult.seed,sep="")
rab$adult.seed2[rab$adult.seed2=="SA"] = "AS"

par(mfrow=c(2,1),mar=c(4,4,1,1))
hist(rab$rab,col="grey",breaks=200,xlab="rab",main = filename[i],xlim=c(0,1))
segments(.8,0,0.8,500,"red")
hist(rab$rab,col="grey",breaks=200,xlab="rab",ylim=c(0,100),main = "",xlim=c(0,1))
segments(.8,0,0.8,100,"red")


### list all pairs that have rab > 0.7; will probably use 0.8 but looking at 0.7
#print(filename[i])
#print(table(rab$adult.seed,rab$rab>0.7))
threshold = 0.7
out = rbind(out,data.frame(
      dataset = filename[i],
      adult.seed = rab$adult.seed2[rab$rab>threshold],
      rab=rab$rab[rab$rab>threshold],
      indA=rab$envA$ind[rab$rab>threshold],
      indB=rab$envB$ind[rab$rab>threshold]))
}
dev.off()

library(reshape)
tmp = cast(melt(out),indA+indB+adult.seed~dataset)
write.table(tmp,"output/clones_r=70perc.txt",quote=F,sep="\t",row.names = F)

### apply 2 criteria
# 1) use only max rab > 0.8
# 2) has to have rab > 0.7 in all four SNP datasets

rabNAs = 4-rowSums(is.na(tmp[,-c(1:3)]))
tmp[is.na(tmp)] = 0
rabMaxs = rowMaxs(tmp[,-c(1:3)])
tmp2 = tmp[rabMaxs > 0.8 & rabNAs == 4,]
write.table(tmp2,"output/clones_r=80perc.txt",quote=F,sep="\t",row.names = F)
print(table(tmp2$adult.seed))
