### Make genotype calls

rm(list=ls())
library(hierfstat)
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


#### only include loci that had NAs in <20% of individuals
prop.ind <- colSums(is.na(geno))/nrow(geno)
geno <- geno[,prop.ind<=0.2]
#> dim(geno)
#[1]   393 8684

write.table(geno,"data/genotypeCalls.393ind.8684loci.geno",sep="\t",quote = F)
