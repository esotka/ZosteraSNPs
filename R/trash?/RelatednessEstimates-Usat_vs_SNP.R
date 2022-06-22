### do usat and SNP based estimates of genetic distance match?

### enter .csv file into genind.
library(poppr)
rm(list=ls())
dat <- read.csv("data/DC_grids_final0521_microsats.csv")
rownames(dat) <- dat$sample; dat <- dat[,-1]
locusname <- colnames(dat)[seq(1,19,2)]
allele1 <- seq(1,19,2); allele2 <- seq(2,20,2)

### metadata genetics
meta <- data.frame(ind=rownames(dat),
                   site=substr(rownames(dat),1,2),
                   site.nice="cur",
                   depth=tolower(substr(rownames(dat),4,4)),
                   grid=substr(rownames(dat),5,5),
                   quad=substr(rownames(dat),7,8))
meta$quad <- ifelse(nchar(meta$quad)==1,paste("0",meta$quad,sep=""),meta$quad)
meta$quad.name <- paste(meta$site.nice,".",meta$depth,meta$grid,"_",meta$quad,sep="")


out <- c()
for (i in 1:length(locusname))
{
  tmp1 <- as.character(dat[,allele1[i]])
  tmp1 <- ifelse(nchar(tmp1)==3,tmp1,paste("0",tmp1,sep=""))
  tmp2 <- as.character(dat[,allele2[i]])
  tmp2 <- ifelse(nchar(tmp2)==3,tmp2,paste("0",tmp2,sep=""))
  out <- cbind(out,paste(tmp1,tmp2,sep=""))
}
rownames(out) <- rownames(dat)
usat <- df2genind(out,ncode=3)

#### get geographic distances for these sites ####

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
m2 <- m[rownames(m)%in%meta$quad.name,colnames(m)%in%meta$quad.name]
m3 <- m2[lower.tri(m2)]

#### Put geno into the same order as geo dist ###
#### calculate relatedness ###

usat2 <- usat[match(rownames(m2),meta$quad.name)]
nei <- nei.dist(usat)
edw <- edwards.dist(usat)
rey <- reynolds.dist(usat)
### make labels
tmp <- matrix(rownames(m2),89,89)
tmp1 <- tmp[lower.tri(tmp)]
tmp2 <- t(tmp)
tmp2 <- tmp2[lower.tri(tmp2)]

#### analysis ###
out <- data.frame(nei.usat=as.vector(nei),
                  edw.usat=as.vector(edw),
                  rey.usat=as.vector(rey),
                  m=m3,ind1=tmp1,ind2=tmp2)


### SNPS ####

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

# subset to the ones that were microsat genotyped
env2 <- env[match(rownames(m2),env$quad.name),]
geno.genind.adult2 <- geno.genind.adult[match(rownames(m2),env$quad.name)]
geno.genind.adult2@ploidy <- as.integer(rep(2,89))
nei <- nei.dist(geno.genind.adult2)
edw <- edwards.dist(geno.genind.adult2)
rey <- reynolds.dist(geno.genind.adult2)

#### get geographic distances for these sites ####

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
m2 <- m[rownames(m)%in%env2$quad.name,colnames(m)%in%env2$quad.name]
m3 <- m2[lower.tri(m2)]


### make labels
tmp <- matrix(rownames(geno.genind.adult2$tab),89,89)
tmp1 <- tmp[lower.tri(tmp)]
tmp2 <- t(tmp)
tmp2 <- tmp2[lower.tri(tmp2)]

#### analysis ###
out.snp <- data.frame(nei.snp=as.vector(nei),
                  edw.snp=as.vector(edw),
                  rey.snp=as.vector(rey))

out.all <- cbind(out,out.snp)

print(xyplot(nei.snp~nei.usat,out.all,type=c("p","r")))



