### glm outlier analysis

rm(list=ls())
library(vcfR)
library(lattice)
dat <- read.vcfR(file="data/SNPs/zos_7chr.214adults.90perc.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])
geno <- t(geno) # row = ind; col = loci
rownames(geno) = colnames(gt2)
colnames(geno) = rownames(gt2)
locNames = colnames(geno)

############# ENV DATA ##############
inds <- rownames(geno)
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
env$site.depth <- paste(env$site.nice,env$SvD,sep="_")

table(env$site.nice,env$SvD)
#         D  S
#  Curlew 39 27
#  Lynch  27 35
#  Niles  23 25
#  West   18 20

p_out = c()
for(i in 1:ncol(geno))
{
# a) remove NAs per locus
tmp = geno[,i]; #names(tmp) = env$SvD
tmp = tmp[complete.cases(tmp)]
depth = substr(names(tmp),2,2)
site = substr(names(tmp),1,1)
# b) glm
p_out = rbind(p_out,data.frame(loc = colnames(geno)[i],site = anova(glm(tmp~site+depth),test="LRT")$"Pr(>Chi)"[2],depth = anova(glm(tmp~site+depth),test="LRT")$"Pr(>Chi)"[3]))
}
outstats = p_out[complete.cases(p_out),]
pdf("output/glm outlier analysis.pdf")
print(histogram(~(-log10(outstats$depth))),breaks=200))
chr <- as.numeric(substr(outstats$loc,4,5))
chr[is.na(chr)] <- 7
cols <- c("red","black","red","black","red","black","red")[chr]
plot(-log10(outstats$depth),ylab="log10 (p-value)",col=cols,pch=20); segments(-10,5.12,10^5,5.12,col="red")
# COrrection: log10(0.05/nrow(outstats))
# [1] -5.120376
outlier = geno[,colnames = outstats$loc[log10(outstats$depth)<(log10(0.05/nrow(outstats)))]]
for (j in 1:ncol(outlier))
{
print(lattice::barchart(prop.table(table(outlier[,j],substr(rownames(outlier),2,2)),2),
stack=F,auto.key = list(space = "right"),main=colnames(outlier)[j]))
a = outlier[,j]
a = a[complete.cases(a)]
depth = substr(names(a),2,2)
site = substr(names(a),1,1)
print(tapply(a,depth,sum)/(2*tapply(a,depth,length)))
}
dev.off()

