### DAPC-based assignment. 
# 1) make groups based on deep vs shallow at each site (8 total)
# 2) assign seeds to these groups

rm(list=ls())
library(adegenet)
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

### adults ###
noClones <- readLines("data/ind245adults_noClones_clean") # adults only 
adult <- out2[rownames(out2)%in%noClones,]
#> dim(geno)
#[1]   245 8684
############# ENV DATA ##############
inds2 <- rownames(adult)
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


dapc2 <- dapc(adult,grp=env$site.nice.depth,n.pca=20,n.da=100)
cols <- c("blue","blue","purple","purple","black","black","red","red")

pdf("output/SeedAssignment-DAPC.pdf")
scatter(dapc2,col=transp(cols,0.7),pch=c(19,15),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=c(1,22),
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=c(1,22),
       cex=3, lwd=2, col=cols)

### seed data ###

AdultSeed <- substr(rownames(out2),7,7)
seed <- geno[AdultSeed=="S",]

############# ENV DATA ##############
inds2 <- rownames(seed)
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

seedPredict <- predict(dapc2,seed)
### no Lynch deep
cols.seed <- c("blue","blue","purple","black","black","red","red")
pch.seed <- c(1,22,22,1,22,1,22)
points(seedPredict$ind.scores,col=cols.seed[env$site.nice.depth],pch=pch.seed[env$site.nice.depth],cex=2)
#par(mar=c(2,2,2,2))
#compoplot(dapc2,new.pred=seedPredict)
out <- as.matrix(table(seedPredict$assign,env$site.nice.depth))
out.perc <- out/table(env$site.nice.depth)
col.l <- colorRampPalette(brewer.pal(9, "Greys"))
f <- levelplot(out,cex=0.2,col.regions=col.l,xlab="Adult meadows",ylab="Dispersed seeds")
print(f)
dev.off()

write.table(out,"output/SeedAssignment-DAPC.txt",sep="\t",quote=F)
