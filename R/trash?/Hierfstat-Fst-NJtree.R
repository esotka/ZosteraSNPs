## do hierfstats FST analysis and Neighbor joining tree with support


rm(list=ls())
library(adegenet)
library(hierfstat)
library(ape)
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

#### pairwise Fst- adults only ###
geno2.all <- data.frame(grp=paste(env$site.nice,env$SvD,sep="_"),geno)
geno2.adults <- geno2.all[env$adult.seed=="A",]
fst <- pairwise.WCfst(geno2.adults) # this takes a long time.
write.csv(fst,"output/fst-adults-SvD.csv")

# make tree
fst.matrix <- read.csv("output/fst-adults-SvD.csv")[,-1]
rownames(fst.matrix) <- colnames(fst.matrix)
tr <- nj(as.matrix(fst.matrix))

# split SNPS into 100 groups (~200 SNPs)
# this takes awhile

geno3.adults <- geno2.adults[,-1] # genotypes only 
snp.breaks <- hist(1:dim(geno3.adults)[2],breaks=100,plot = F)$breaks
tr.100reps <- list()
  
for (i in 1:98)
{
  print(i)
  tmp <- data.frame(grp=geno2.adults[,1],geno3.adults[,(snp.breaks[i]+1):snp.breaks[i+1]])
  tmpfst <- pairwise.WCfst(tmp) 
  tr.100reps[[i]] <- nj(tmpfst)
}
save(x = tr.100reps,file = "output/trees100reps.RData")

load("output/trees100reps.RData")
# 96 groups of 200 SNPs each
print(prop.part(tr.100reps[-97])) #use this output to make statistical support of branches

pdf("output/fst-adults-SvD-NJtree.pdf")
plot(tr,"unrooted")
add.scale.bar(x=0.004,y=0.005)
dev.off()
