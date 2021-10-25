### convert vcf of genotype calls from bcftools view (CC, CT, etc..) into single letter IUPAC codes and then make a fasta file
rm(list=ls())
library(vcfR) # read.vcfR()
library(seqinr) # bma() Computing an IUPAC nucleotide symbol
dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=T)
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

inds <- readLines("data/ind393_clean")
loci.to.use <- read.table("data/loci19433")
loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")

### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use$V3,]

library(parallel)

## by locus.

out <- mclapply(1:nrow(gt2), function(i)
  {
tmp1 <- paste(substr(gt2[i,],1,1),substr(gt2[i,],3,3),sep="")
tmp1 <- ifelse(nchar(tmp1)==1,"N",tmp1)
tmp2 <- mclapply(as.list(tmp1),FUN=function(j) bma(s2c(j)))
})

out2 <- as.data.frame(matrix(unlist(out),nrow=393)) # 393 individuals
out2.noNAs <- out2[,colSums(is.na(out2))==0] ### get rid of all NAs

write.table(x = out2.noNAs,"data/fasta/allLoci.txt",quote=F,row.names = F,col.names = F,sep = "")

# now make a fasta file

for (i in 1:length(inds))
{
  tmp <- read.table("data/fasta/allLoci.txt",skip=i-1,nrows=1)
  if (i == 1)
  {
    write.table(x = paste(">",inds[i],"\n",tmp,"\n",sep=""),"data/fasta/393ind.fas",quote=F,row.names = F,col.names = F)
  }
  else
  {
    write.table(x = paste(">",inds[i],"\n",tmp,"\n",sep=""),"data/fasta/393ind.fas",quote=F,row.names = F,col.names = F,append = T)
  }
}


