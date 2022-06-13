### PCA on genotype calls from bcftools

rm(list=ls())
library(adegenet)
library(vcfR)

dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
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

## make an noClones + seeds dataset
inds <- readLines("data/ind393_clean")
seeds <- inds[substr(inds,7,7)=="S"]
noClones <- readLines("data/ind245adults_noClones_clean")
geno <- geno[rownames(geno)%in%c(noClones,seeds),]
#> dim(geno)
#[1]   344 19432

#### only include loci that had NAs in <50% of individuals
prop.ind <- colSums(is.na(geno))/nrow(geno)
geno <- geno[,prop.ind<=0.5]
#> dim(geno)
#[1]   245 12951

### replace NAs with the most common genotype per locus
library(parallel)

out <- mclapply(1:ncol(geno), function(i)
{
  tmp <- geno[,i]
  common.gt <- names(table(tmp))[table(tmp)==max(table(tmp))]
  tmp[is.na(tmp)] <- common.gt
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(geno)))


pdf('output/PCA_genotypeCalls.pdf')
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
env$SvD <- factor(env$SvD)
cols <- c("blue","purple","black","red")

# all
pca1 <- dudi.pca(out2,center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca1$li,env$site.nice,col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults + seeds",line=2)
# adults
pca2 <- dudi.pca(out2[env$adult.seed=="A",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca2$li,env$site.nice[env$adult.seed=="A"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults only",line=2)

# seeds
pca3 <- dudi.pca(out2[env$adult.seed=="S",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca3$li,env$site.nice[env$adult.seed=="S"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("seeds",line=1)

## adult SvD
for (i in 1:4)
{
pca4 <- dudi.pca(out2[env$site.nice==levels(env$site.nice)[i]  & env$adult.seed=="A",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca4$li,env$SvD[env$site.nice==levels(env$site.nice)[i]  & env$adult.seed=="A"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext(levels(env$site.nice)[i],line=2)
}

dev.off()


