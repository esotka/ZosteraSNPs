### correlogram - genetic data vs distance
### https://leandroroser.github.io/EcoGenetics-Tutorial/
### also borrowed from https://bookdown.org/hhwagner1/LandGenCourse_book/WE-5.html 
# library(devtools)
# install_github("leandroroser/EcoGenetics-devel")
rm(list=ls())
library(EcoGenetics)
library(hierfstat)
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


############# ENV DATA ##############
inds <- readLines("data/ind393_clean")
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.depth.lifeHistory <- paste(env$site.nice,env$SvD,env$adult.seed,sep="_")
env$position <- as.numeric(substr(env$ind,5,6))
env$site.depth.position <- paste(env$site,env$SvD,env$position,sep="_")

### geographic coordinates
geo <- read.csv("output/IndivLatLon.csv")
rownames(geo) <- geo[,1]
geo$pop[geo$pop=="cur"] <- "D"
geo$pop[geo$pop=="nil"] <- "N"
geo$pop[geo$pop=="wes"] <- "W"
geo$pop[geo$pop=="lyn"] <- "L"
geo$SvD <- toupper(substr(geo$grid,5,5))
geo$site.depth.position <- paste(geo$pop,geo$SvD,geo$ind,sep="_")

## adults only ##
env.ad <- env[env$adult.seed=="A",]
geno.ad <- geno[env$adult.seed=="A",]
geo.ad <- geo[match(env.ad$site.depth.position,geo$site.depth.position),]

#### correlogram ###

site.depth.unique <- unique(env.ad$site.depth.lifeHistory)

for (i in 1:length(site.depth.unique))
{
  env.ad.tmp <- env.ad[env.ad$site.depth.lifeHistory==site.depth.unique[i],]
  geno.ad.tmp <- geno.ad[env.ad$site.depth.lifeHistory==site.depth.unique[i],]
  geo.ad.tmp <- geo.ad[env.ad$site.depth.lifeHistory==site.depth.unique[i],c("lon","lat")]
  geo.ad.tmp$lon <- -(geo.ad.tmp$lon)
  rownames(geo.ad.tmp) <- env.ad.tmp$ind
  nas <- colSums(is.na(geno.ad.tmp)) # number of NAs per locus
  #95% threshold
  loci.to.use <- nas<=(.05*(nrow(geno.ad.tmp)))
  geno.ad.tmp2 <- geno.ad.tmp[,loci.to.use]
  eco.tmp <- ecogen(XY=geo.ad.tmp) # make a ecogen object
  ecoslot.G(eco.tmp,order.G = TRUE) <- geno.ad.tmp2
  moran <- eco.correlog(Z = eco.tmp[["A"]], XY = eco.tmp[["XY"]],latlon=T,method="I",size=100,int=3) # this worked 
}
size=100,int=2,
# create eco frame
eco.temp <- ecogen(G = geno)

eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, E = environment, S = structure, order.G = TRUE)


moran.obs <- lapply(moran@OUT,"[[",2)
xbar <- c()
for (mm in 1:8){xbar <- c(xbar,mean(unlist(lapply(moran.obs,"[[",mm)),na.rm=T))}
