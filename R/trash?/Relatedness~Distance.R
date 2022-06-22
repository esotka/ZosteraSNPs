### relatedness with distance - Zostera
rm(list=ls())
library(adegenet)
library(poppr)
geno <- read.table('data/zos.393ind.HWE.99.calls.forGenind',sep=" ")
geno <- t(geno)
colnames(geno) <- readLines("data/loci19433")
rownames(geno) <- readLines("data/ind393_clean")

geno[geno==0] <- "11" 
geno[geno=="1"] <- "12"
geno[geno=="2"] <- "22"

geno.genind <- df2genind(geno,ncode=1,NA.char=NA)
names(geno.genind@all.names) <- colnames(geno)

### adults only ###
inds <- readLines("data/ind393_clean")
geno.genind.adult <- geno.genind[substr(inds,7,7)=="A"]

### meta of adults####
inds.adult <- inds[substr(inds,7,7)=="A"]
env <- data.frame(ind=substr(inds.adult,1,7),
                  adult.seed=substr(inds.adult,7,7),
                  site=substr(inds.adult,1,1),
                  SvD=tolower(substr(inds.adult,2,2)),
                  grid=substr(inds.adult,3,3),
                  ind2=substr(inds.adult,5,6))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "cur"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "nil"
env$site.nice[env$site.nice=="W"] <- "wes"
env$site.nice[env$site.nice=="L"] <- "lyn"
env$site.nice <- factor(env$site.nice)
env$quad.name <- paste(env$site.nice,".",env$SvD,env$grid,"_",env$ind2,sep="")


#### get geographic distances for these sites ####

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
m2 <- m[rownames(m)%in%env$quad.name,colnames(m)%in%env$quad.name]
m3 <- m2[lower.tri(m2)]

#### Put geno into the same order as geo dist ###
#### calculate relatedness ###

geno.genind.adult2 <- geno.genind.adult[match(rownames(m2),env$quad.name)]
nei <- nei.dist(geno.genind.adult2)
edw <- edwards.dist(geno.genind.adult2)
rey <- reynolds.dist(geno.genind.adult2)
#rog <- rogers.dist(geno.genind.adult2) # takes a long time

### make labels
tmp <- matrix(rownames(m2),294,294)
tmp1 <- tmp[lower.tri(tmp)]
tmp2 <- t(tmp)
tmp2 <- tmp2[lower.tri(tmp2)]

#### analysis ###
out <- data.frame(nei=as.vector(nei),
                  edw=as.vector(edw),
                  rey=as.vector(rey),
                  m=m3,ind1=tmp1,ind2=tmp2)
out$site <- ifelse(substr(out$ind1,1,3)==substr(out$ind2,1,3),substr(out$ind1,1,3),"NA")
out$depth <- ifelse(substr(out$ind1,5,5)==substr(out$ind2,5,5),substr(out$ind1,5,5),"NA") 
out <- out[!(out$site=="NA" | out$depth=="NA"),]

library(lattice)
pdf("output/Relatedness~Distance.pdf")
dist.names <- c("nei","edw","rey")
for(i in 1:3)
{
print(xyplot(out[,dist.names[i]]~m,type=c("p","r"),data=out,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
  print(xyplot(out[,dist.names[i]]~m | site+depth,type=c("p","r"),data=out,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
}
dev.off()

### stats ####
out$site.depth <- paste(out$site,out$depth,sep="_")
site.depth.unique <- unique(out$site.depth)
out.stats <- c()
for (j in 1:length(site.depth.unique))
{
for(i in 1:3)
{
tmp <- data.frame(x=out$m[out$site.depth==site.depth.unique[j]]
                  ,y=out[out$site.depth==site.depth.unique[j],dist.names[i]])
tmp.stats <- cor.test(tmp$x,tmp$y,method="kendall")
out.stats <- rbind(out.stats,
                   data.frame(site.depth=site.depth.unique[j],
                              dist.names=dist.names[i],
                              n=dim(tmp)[1],
                              stat=tmp.stats$statistic,
                              p=round(tmp.stats$p.value,3),
                              signif=ifelse(tmp.stats$p.value<0.05,"*","")))
}
}
print(out.stats)
  
