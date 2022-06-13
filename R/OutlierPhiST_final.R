### outlier PhiST

rm(list=ls())
library(pegas)
library(ade4)
library(hierfstat)

geno <- read.delim("data/genotypeCalls.393ind.8684loci.geno")

noClones <- readLines("data/ind245adults_noClones_clean") # adults only 
geno <- geno[rownames(geno)%in%noClones,]
#> dim(geno)
#[1]   245 8684

############# ENV DATA ##############
inds2 <- rownames(geno)
env <- data.frame(ind=inds2,
                  adult.seed=substr(inds2,7,7),
                  site=substr(inds2,1,1),
                  SvD=substr(inds2,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
env$site.nice.depth <- factor(paste(env$site.nice,env$SvD,sep="_"))

geno2.all <- data.frame(grp=env$site.nice.depth,geno)

for (i in 1:length(levels(env$site.nice)))
{
  tmp <- geno2.all[env$site.nice==levels(env$site.nice)[i],]
  tmp.fst <- pairwise.WCfst(tmp[,c(1,2)])[2]
  
  for(j in 3:dim(tmp)[2])
  {  
    if (length(table(tmp[,j]))==1)
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
loci <- readLines("data/loci8684")
chr <- unlist(lapply(strsplit(loci,"_"),"[[",1))
out <- data.frame(site=c(rep("Curlew",8684),rep("Lynch",8684),rep("Niles",8684),rep("West",8684)),
                  loci=rep(readLines("data/loci8684"),4),
                  chr=rep(chr,4),
                  phist)

phist.sub <- c()
for(k in 1:4)
{
  tmp <- out[out$site==unique(out$site)[k],]
  tmp <- tmp[order(tmp$phist,decreasing = T),]
  phist.sub <- rbind(phist.sub,tmp[1:87,]) # top 1%
}

### result ==> there are no SNPs that are outliers (1% of highest PhiSt) in all four locations. There are no SNPs that are outliers in 3 locations. There are two SNPs that are outliers in 2 of 4 locations

table(table(phist.sub$loci))
#1   2 
#330   9 

### outlier Fst plot

pdf("output/outlierPhiSt.pdf")
par(mfrow=c(4,1),mar=c(1,2,2,0))
for(k in 1:4)
{
  tmp <- out[out$site==unique(out$site)[k],]
  plot(tmp$phist,cex=.1,pch=20,ylim=c(0,.5),xaxt="n")
  mtext(unique(out$site)[k],line=-2)
}
dev.off()


