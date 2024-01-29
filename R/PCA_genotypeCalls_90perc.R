### PCA on genotype calls from bcftools

rm(list=ls())
library(adegenet)
library(vcfR)

dat = read.vcfR("data/SNPs/zos_7chr.214adults.99seeds.90perc.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])
geno <- t(geno) # row = ind; col = loci

rownames(geno) = colnames(gt2)
colnames(geno) = rownames(gt2)
inds = rownames(geno)
#> dim(geno)
#[1]   344 6611

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


pdf("output/PCA_genotypeCalls_90perc.pdf")
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
env$SvD.adult.seed = factor(paste(env$SvD,env$adult.seed,sep=""))
env$site.nice.depth = factor(paste(env$site.nice,env$SvD))
cols <- c(rep("red",2),rep("darkgrey",2),rep("black",2),rep("forestgreen",2))
pchs = c(19,3,21,3)

#adults + seeds 
pca1 <- dudi.pca(out2,center=TRUE,scale=FALSE,scannf=FALSE,nf=2)
s.class(pca1$li,env$site.nice.depth,
      col=transp(cols,.6),grid=F,cstar = 0,cpoint=2,
      pch=pchs[env$SvD.adult.seed])

legend(x=10,y=-7,c("Seed","Deep Adult","Shallow Adult"),pch=c(3,19,21))
dev.off()




