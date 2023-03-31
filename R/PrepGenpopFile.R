## data for NeEstimator in Genpop

rm(list=ls())
library(adegenet)
library(RColorBrewer)
library(lattice)
geno <- read.delim("data/genotypeCalls.393ind.8684loci.geno")
### replace NAs with most common genotype
library(parallel)

out <- mclapply(1:ncol(geno), function(i)
{
  tmp <- geno[,i]
  common.gt <- names(table(tmp))[table(tmp)==max(table(tmp))]
  tmp[is.na(tmp)] <- common.gt
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(geno))) # rows = 393 ind cols = loci
rownames(out2) <- rownames(geno)
colnames(out2) <- colnames(geno)
### seed data ###

AdultSeed <- substr(rownames(out2),7,7)
seed <- out2[AdultSeed=="S",]

seed[seed==0] <- "11"
seed[seed==1] <- "12"
seed[seed==2] <- "22"

write.table(seed,"output/seed_genpop.txt",quote=F,rownames=F,colnames=F)
writeLines(rownames(seed),"output/seed_names")
