### Are there clones?

rm(list=ls())
library(adegenet)
library(poppr)
library(adegenet)
library(vcfR)

dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- "11"
gt2[gt2=="0/1"] <- "12"
gt2[gt2=="1/1"] <- "22"
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

loci.to.use <- read.table("data/loci19433")
loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")
### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use$V3,]
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])
geno <- t(geno) # row = ind; col = loci

rownames(geno) <- readLines("data/ind393_clean")
colnames(geno) <- rownames(gt2)


### meta of adults####
inds <- readLines("data/ind393_clean")
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
env$site.nice.depth <- paste(env$site.nice,env$SvD,sep="_")
env$site.nice.depth <- factor(env$site.nice.depth)


## adults only
geno <- geno[rownames(geno)%in%env$ind,]

### number of genotypes across populations
x <- new("genlight",geno,ploidy=2)
x@pop <- env$site.nice.depth
x2 <- as.snpclone(x)
x2@pop <- env$site.nice.depth

mlgID <- list()
threshold = 0.01
pdf("output/MLG_analysis-threshold=0.01.pdf",width=10,height=4)
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1), c(2,1,0,0))
for (i in 1:8)
{
x2.sub <- x2[x2@pop==levels(env$site.nice.depth)[i]]
bwd <- bitwise.dist(x2.sub,missing_match=T)
hist(bwd,col="grey",breaks=100,xlim=c(0,0.1),main=levels(env$site.nice.depth)[i])
segments(threshold,-1,threshold,100,col="red")
# Visual inspection suggests threshold of 0.01
mlgID[[i]] <- mlg.filter(x2.sub,threshold = threshold,distance=bwd)
}


names(mlgID) <- levels(env$site.nice.depth)

## frequency of unique MLGs
print(num.uniqueMLG <- unlist(lapply(lapply(mlgID,unique),length)))
print(sampleSize <- unlist(lapply(mlgID,length)))
print(freq.uniqueMLG <- num.uniqueMLG/sampleSize)

# Where are the clones? 
# here are geographic positions
latlon <- read.csv("output/IndivLatLon.csv")
tmp <- ifelse(latlon$ind<10,paste("0",latlon$ind,sep=""),latlon$ind)
latlon$ind <- tmp
latlon$quad.name <- paste(latlon$grid,latlon$ind,sep = "_")

for (i in 1:8)
{
cloneID <- mlgID[[i]]
tmpID <- env$ind[env$site.nice.depth==levels(env$site.nice.depth)[i]]
tmpQuadName <- env$quad.name[env$site.nice.depth==levels(env$site.nice.depth)[i]]
tmpLat <- latlon$lat[match(tmpQuadName,latlon$quad.name)]
tmpLon <- latlon$lon[match(tmpQuadName,latlon$quad.name)]
all <- data.frame(cloneID,tmpID,tmpQuadName,tmpLat,tmpLon)
plot(tmpLon,tmpLat,main=levels(env$site.nice.depth)[i])
### colored are clones
cloneID.plot <- names(table(cloneID))[table(cloneID)>1]
print(cloneID.plot)
if(length(cloneID.plot)==1){
  all.sub <- all[all$cloneID==cloneID.plot,]
  for (k in 1:dim(all.sub)[1])
  {
    if (dim(all.sub)[1]==2) {segments(all.sub$tmpLon[1],all.sub$tmpLat[1],all.sub$tmpLon[2],all.sub$tmpLat[2],col="red")}
    else {
    x <- combn(all.sub$tmpLon,2)
    y <- combn(all.sub$tmpLat,2)
    segments(x[1,k],y[1,k],x[2,k],y[2,k],col="red")}}
}
if(length(cloneID.plot)==2)
  {
  all.sub <- all[all$cloneID==cloneID.plot[1],]
  for (k in 1:dim(all.sub)[1])
  {
    if (dim(all.sub)[1]==2) {segments(all.sub$tmpLon[1],all.sub$tmpLat[1],all.sub$tmpLon[2],all.sub$tmpLat[2],col="blue")}
    else {
    x <- combn(all.sub$tmpLon,2)
    y <- combn(all.sub$tmpLat,2)
    segments(x[1,k],y[1,k],x[2,k],y[2,k],col="red")}}
  all.sub <- all[all$cloneID==cloneID.plot[2],]
  for (k in 1:dim(all.sub)[1])
  {
    if (dim(all.sub)[1]==2) {segments(all.sub$tmpLon[1],all.sub$tmpLat[1],all.sub$tmpLon[2],all.sub$tmpLat[2],col="blue")}
        else {
          x <- combn(all.sub$tmpLon,2)
          y <- combn(all.sub$tmpLat,2)
          segments(x[1,k],y[1,k],x[2,k],y[2,k],col="blue")
        }
    }
}

}


dev.off()




