### ngsRelate results
rm(list=ls())
### meta of adults####
inds.adult <- readLines("data/ind294adults_clean")
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


#### get geographic distances for these sites ####

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
m2 <- m[rownames(m)%in%env$quad.name,colnames(m)%in%env$quad.name]
m3 <- m2[lower.tri(m2)]

### ngsRelate  "rab"
### commands from ngsRelate
### see https://github.com/ANGSD/NgsRelate
### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
#angsd -b DD_43_bamfiles -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -rf zos.393ind.HWE.99.forANGSD.loci19433.txt -out DD
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
#zcat DD.mafs.gz | cut -f5 |sed 1d > DD_freq
### run NgsRelate The output should be a file called newres that contains the output for all pairs
#./ngsRelate -g DD.glf.gz -n 43 -f DD_freq -O DD_newres


a <- c("DD",
       "DS",
       "LD",
       "LS",
       "ND",
       "NS",
       "WD",
       "WS")

rab <- c()
for (i in 1:length(a))
{
  tmp <- read.delim(file=paste("data/ngsRelate/",a[i],"_newres",sep=""))[,c("a","b","rab")]
  tmpind <- readLines(paste("data/ngsRelate/",a[i],sep=""))
  tmp.quad.name <- env$quad.name[match(tmpind,env$ind)]
  tmp.key <- data.frame(num=0:(length(tmpind)-1),tmp.quad.name)
  tmp$a <- tmp.key$tmp.quad.name[match(tmp$a,tmp.key$num)]
  tmp$b <- tmp.key$tmp.quad.name[match(tmp$b,tmp.key$num)]
  rab <- rbind(rab,tmp)
}

rab$site <- substr(rab$a,1,3)
rab$depth <- substr(rab$a,5,5)
rab$site_depth <- paste(rab$site,rab$depth,sep="_")
library(lattice)
print(histogram(~rab | site + depth,data=rab,col="grey",breaks=40))
print(histogram(~log(rab+0.000001) | site + depth,data=rab,col="grey",breaks=40))

## get geographic distance
m.all <- c()
for(j in 1:dim(rab)[1])
{m.all <- c(m.all,m2[rownames(m2)==rab$a[j],colnames(m2)==rab$b[j]])}
rab$m <- m.all
pdf("output/Relatedness~Distance-ngsRelate.pdf")
print(xyplot(rab~m | site+depth,data=rab,cex=0.5,pch=20,lwd=4,col="grey",type=c("p","r")))

dev.off()

### stats - Kendall's Tau
site.depth.unique <- unique(rab$site_depth)
out.stats <- c()
for (j in 1:length(site.depth.unique))
{
    tmp <- data.frame(x=rab$m[rab$site_depth==site.depth.unique[j]]
                      ,y=rab$rab[rab$site_depth==site.depth.unique[j]])
    tmp.stats <- cor.test(tmp$x,tmp$y,method="kendall")
    out.stats <- rbind(out.stats,
                       data.frame(site.depth=site.depth.unique[j],
                                  n=dim(tmp)[1],
                                  stat=tmp.stats$statistic,
                                  p=round(tmp.stats$p.value,3),
                                  signif=ifelse(tmp.stats$p.value<0.05,"*","")))
}
print(out.stats)

### stats - randomization ANCOVA
### code borrowed from Allan Strand (CofC)

### data are zero-inflated.
### Used a glm with Poisson distribution for each site independently
### AER::dispersiontest() indicated this was not overdispersed (p=1.000)
### Interaction between site, depth and geographic distance is significant
### Used pgirmess::PermTest() for permutation test.
#Patrick Giraudoux (2021). pgirmess: Spatial Analysis and Data Mining
#for Field Ecologists. R package version 1.7.0.
#https://CRAN.R-project.org/package=pgirmess

mod <- glm(rab~m*site*depth,data=rab,family="poisson")
library(AER)
dispersiontest(mod)
anova(mod)

#Df Deviance Resid. Df Resid. Dev
#NULL                          5422     664.22
#m             1   82.135      5421     582.08
#site          3   10.492      5418     571.59
#depth         1    2.467      5417     569.12
#m:site        3    4.540      5414     564.58
#m:depth       1    1.849      5413     562.73
#site:depth    3    1.177      5410     561.56
#m:site:depth  3    6.736      5407     554.82

#print(PermTest(mod,B=1000)) ### takes awhile

#Call: 
#  PermTest.glm(obj = mod, B = 1000)

#Based on 1000 replicates
#Simulated p-value:
#  p.value
#m              0.000
#site           0.000
#depth          0.001
#m:site         0.000
#m:depth        0.006
#site:depth     0.179
#m:site:depth   0.000

# Analyzed m*depth at four sites
#for (i in 1:4)
#{
#  tmp <- rab[rab$site==unique(rab$site)[i],]
#  mod.tmp <- glm(rab~m*depth,data=tmp,family="poisson")
#  out <- PermTest(mod.tmp,B=1000)
#  print(unique(rab$site)[i])
#  print(anova(mod.tmp))
#  print(out)
#}

#[1] "cur"
#
#        Df Deviance Resid. Df Resid. Dev
#NULL                     1682     197.18
#m        1  22.9682      1681     174.21
#depth    1   2.9017      1680     171.31
#m:depth  1   7.9286      1679     163.38
#
#Based on 1000 replicates
#Simulated p-value:
#        p.value
#m             0
#depth         0
#m:depth       0
#
#[1] "lyn"
#
#        Df Deviance Resid. Df Resid. Dev
#NULL                     1848     204.02
#m        1  30.8058      1847     173.22
#depth    1   0.4486      1846     172.77
#m:depth  1   0.0238      1845     172.74
#
#Based on 1000 replicates
#Simulated p-value:
#        p.value
#m         0.000
#depth     0.171
#m:depth   0.743
#
#[1] "nil"
#
#Terms added sequentially (first to last)
#
#
#        Df Deviance Resid. Df Resid. Dev
#NULL                     1155     124.62
#m        1   8.7298      1154     115.89
#depth    1   0.0080      1153     115.88
#m:depth  1   0.0002      1152     115.88
#
#Based on 1000 replicates
#Simulated p-value:
#        p.value
#m         0.000
#depth     0.851
#m:depth   0.969
#
#[1] "wes"
#
#        Df Deviance Resid. Df Resid. Dev
#NULL                      734     128.39
#m        1  24.7371       733     103.66
#depth    1   0.0540       732     103.60
#m:depth  1   0.7849       731     102.82
#
#Based on 1000 replicates
#Simulated p-value:
#        p.value
#m         0.000
#depth     0.689
#m:depth   0.122


site.depth.unique <- unique(rab$site_depth)
out.stats <- c()
pdf("output/Relatedness~Distance-ngsRelate-final.pdf",height=6,width=10)
par(mfcol=c(2,4),mai=c(.3,.3,.1,.1), c(2,1,0,0))
for (j in 1:length(site.depth.unique))
{
  tmp <- data.frame(x=rab$m[rab$site_depth==site.depth.unique[j]]
                    ,y=rab$rab[rab$site_depth==site.depth.unique[j]])
  mod.geo <- glm(y~x,data=tmp,family="poisson")
  x.predict <- seq(1,45,0.5)
  rab.predict = predict(mod.geo,list(x=x.predict),type="response")
  plot(tmp$y~tmp$x,xlim=c(0,40),ylim=c(0,1),ylab="",xlab="")#,xaxt="n",yaxt="n")
  lines(x.predict,rab.predict,col="red",lwd=4)
  mtext(site.depth.unique[j],line=-2)
}
dev.off()
