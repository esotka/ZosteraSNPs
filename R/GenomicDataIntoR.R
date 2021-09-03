### take called genotypes from angsd and output a genind object

## Method
rm(list=ls())
library(adegenet)
dat <- read.delim('data/zos.393ind.HWE.99.geno.geno.gz',header=F)
dat <- dat[,-dim(dat)[2]] # last column is NA
dat[dat=="NN"] <- NA
nind <- 393


for (i in 1:nrow(dat)) #row = SNPs
{
tmp <- as.character(dat[i,-c(1:4)])
majAllele <- dat[i,3]
minAllele <- dat[i,4]
geno3 <- data.frame(geno=c(paste(majAllele,majAllele,sep=""),
                           paste(majAllele,minAllele,sep=""),
                           paste(minAllele,minAllele,sep="")),
                    score=c(2,1,0))
tmp2 <- geno3$score[match(tmp,geno3$geno)]
if(i == 1)
{
  write(tmp2,"data/zos.393ind.HWE.99.calls.forGenind",ncolumns = nind)#,sep=" ",row.names = F,col.names = F,quote=F)}
}

else
  {write(tmp2,"data/zos.393ind.HWE.99.calls.forGenind",append=T,ncolumns = nind)#,row.names = F,col.names = F,quote=F)
}}

out <- read.delim("data/zos.393ind.HWE.99.calls.forGenind",header=F,sep=" ")
geno <- t(out) # row = ind; col = SNP
dat$V1 <- gsub(".1","",dat$V1) ### get rid of .1
colnames(geno) <- paste(dat$V1,dat$V2,sep="_")
writeLines(colnames(geno),"data/loci19433")
id <- readLines("data/ind393")
id[id=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"
rownames(geno) <- substr(id,1,7)
write.table(geno,"data/zos.393ind.HWE.99.calls.forGenind2.txt",quote=F)

#geno2 <- read.delim("zos.393ind.HWE.99.calls.forGenind2.txt",sep=" ")
geno.genind <- df2genind(geno,ncode=1,NA.char=NA)
names(geno.genind@all.names) <- colnames(geno)


### check that it worked with DAPC
gpa <- find.clusters(geno.genind,n.clust = 3,n.pca = 2)
dapc1 <- dapc(geno.genind,gpa$grp,n.pca = 2,n.da = 10)
scatter(dapc1)


