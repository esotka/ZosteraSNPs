### create files for running spagedi 
### 1) pull in genotypes from bcftools (vcf). subset to the loci that occur in 95% of samples (for time expediency)
### 2) pull in geographic distances already calculated. 
### 3) make names consistent.
### 4) use EcoGenetics::ecogen2spagedi() to create 

rm(list=ls())
#library(EcoGenetics)
#library(hierfstat)
library(vcfR)
dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 11
gt2[gt2=="0/1"] <- 12
gt2[gt2=="1/1"] <- 22
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

## nonclonal adults only ##
inds <- readLines("data/ind393_clean")
noClones <- inds%in%readLines("data/ind245adults_noClones_clean")
geno <- geno[noClones,]

# geographic distance in column form
inds.adult <- readLines("data/ind393_clean")
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

### non clonal adults
env <- env[env$ind%in%readLines("data/ind245adults_noClones_clean"),]

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
library(reshape2)
m[upper.tri(m)] <- NA
diag(m) <- NA
m2 <- as.matrix(m)
m.forSpaghedi <- melt(m2)
m.forSpaghedi <- m.forSpaghedi[complete.cases(m.forSpaghedi$value),]
m.forSpaghedi <- m.forSpaghedi[m.forSpaghedi$Var1%in%env$quad.name & m.forSpaghedi$Var2%in%env$quad.name,]
m.forSpaghedi$Var1 <- env$ind[match(m.forSpaghedi$Var1,env$quad.name)]
m.forSpaghedi$Var2 <- env$ind[match(m.forSpaghedi$Var2,env$quad.name)]

### save files

site.depth.unique <- unique(env$site.nice.depth)

for (i in 1:length(site.depth.unique))
{
  
  env.tmp <- env[env$site.nice.depth==site.depth.unique[i],]
  geno.tmp <- geno[env$site.nice.depth==site.depth.unique[i],]
  nas <- colSums(is.na(geno.tmp)) # number of NAs per locus
  #75% threshold
  loci.to.use <- nas<=(.75*(nrow(geno.tmp)))
  geno.tmp2 <- geno.tmp[,loci.to.use]
  geno.tmp2[is.na(geno.tmp2)] <- "00"
  print(site.depth.unique[i])
  print(dim(geno.tmp2))
  geno.tmp2 <- data.frame(ind=rownames(geno.tmp2),pop=1,geno.tmp2)
  write.table(geno.tmp2,paste("~/Desktop/",site.depth.unique[i],"_input.txt",sep=""),quote=F,row.names=F,col.names = T,sep="\t")
  m.forSpaghedi.tmp <- m.forSpaghedi[m.forSpaghedi$Var1%in%rownames(geno.tmp) & m.forSpaghedi$Var2%in%rownames(geno.tmp),]
  write.table(m.forSpaghedi.tmp,paste("~/Desktop/",site.depth.unique[i],"_geo.txt",sep=""),quote=F,row.names=F,col.names = F,sep="\t")
  
}

## Manually add the other stuff to genotype file
#// #ind #cat #coord #loci #dig/loc ploidy
# 40	0	0	10032	2	2  # make sure these are tab-deliminated
# -8
# ind pop <names of loci>
# ...
# END

### Manually add other stuff to distance file
# C820
# ...

#END

## run spagedi with a command file "spagedi < cmds.txt"
# Curlew_D_A_input.txt
# Curlew_D_A_out.txt
#
# 1
# 13
# Curlew_D_A_geo.txt
# 2
# 100
# 100




# test with 75% threshold#

dat <- read.csv("~/Desktop/spagedi_analysis/Curlew_S_A_out.csv")
pdf("~/Desktop/spagedi_analysis/Curlew_S_Ad.spagedi.pdf")
plot(x=dat$Mean.distance,y=dat$Obs.val,xlab="m",ylab="kinship",type="b",ylim=c(-.1,.2))
points(x=dat$Mean.distance,y=dat$X95.CI.inf,col="red",lty="dashed",type="l")
points(x=dat$Mean.distance,y=dat$X95.CI.sup,col="red",lty="dashed",type="l")
text(dat$Mean.distance,y=.2,ifelse(dat$P.2.sided.test..H1..obs..exp.<0.05,"**","ns"))
segments(-1,0,50,0)
dev.off()







