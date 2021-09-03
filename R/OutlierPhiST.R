## do hierfstats FST analysis and Neighbor joining tree with support


rm(list=ls())
library(adegenet)
library(hierfstat)
geno <- read.delim("data/zos.393ind.HWE.99.calls.forGenind",sep=" ",header=F)
geno <- t(geno) # row = ind; col = loci
inds <- readLines("data/ind393")
inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"
inds <- substr(inds,1,7)
#names(geno.genind@all.names) <- inds
############# ENV DATA ##############
env <- data.frame(ind=inds,
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)


### pairwise Fst - by locus - between Deep and Shallow ecotypes
## this takes awhile

geno2.all <- data.frame(grp=paste(env$site.nice,env$SvD,sep="_"),geno)

for (i in 1:length(levels(env$site.nice)))
{
tmp <- geno2.all[env$site.nice==levels(env$site.nice)[i],]
tmp.fst <- pairwise.WCfst(tmp[,c(1,2)])[2]

for(j in 3:dim(tmp)[2])
{  
  if (length(unique(tmp[,j]))==1)
  {tmp.fst <- c(tmp.fst,0)}
  else {
    tmp.fst <- c(tmp.fst,pairwise.WCfst(tmp[,c(1,j)])[2])
  }
}
  write.table(tmp.fst,file = paste("output/PhiST-",levels(env$site.nice)[i]),row.names=F,col.names=F)
}

### 99% outliers

phist <- as.numeric(c(readLines("output/PhiST- Curlew"),readLines("output/PhiST- Lynch"),readLines("output/PhiST- Niles"),readLines("output/PhiST- West")))
phist <- ifelse(phist<=0,0,phist)
loci <- readLines("data/loci19433")
chr <- unlist(lapply(strsplit(loci,"_"),"[[",1))
out <- data.frame(site=c(rep("Curlew",19433),rep("Lynch",19433),rep("Niles",19433),rep("West",19433)),
                  loci=rep(readLines("data/loci19433"),4),
                  chr=rep(chr,4),
                  phist)

phist.sub <- c()
for(k in 1:4)
{
  tmp <- out[out$site==unique(out$site)[k],]
  tmp <- tmp[order(tmp$phist,decreasing = T),]
  phist.sub <- rbind(phist.sub,tmp[1:194,])
}

outlier3sites <- names(table(phist.sub$loci))[table(phist.sub$loci)==3]


### result ==> there are three loci that are outliers (1% of highest PhiSt) in at least three locations. 
 print(phist.sub[phist.sub$loci%in%outlier3sites,])

       site            loci      chr      phist
3727  Curlew LFYR0006_336140 LFYR0006 0.08288411
11640 Curlew    LFYR05_72696   LFYR05 0.04876252
11641 Curlew    LFYR05_72706   LFYR05 0.04876252
31073  Lynch    LFYR05_72696   LFYR05 0.06493195
31074  Lynch    LFYR05_72706   LFYR05 0.06493195
50506  Niles    LFYR05_72696   LFYR05 0.03317238
50507  Niles    LFYR05_72706   LFYR05 0.03317238
42593  Niles LFYR0006_336140 LFYR0006 0.02371256
62026   West LFYR0006_336140 LFYR0006 0.05435535

> tmp <- phist.sub[phist.sub$loci%in%outlier3sites,]
> table(tmp$site,tmp$loci)

           LFYR0006_336140 LFYR05_72696 LFYR05_72706
Curlew               1            1            1
Lynch                0            1            1
Niles                1            1            1
West                 1            0            0


