### convert beagle formatted file of genotype calls (CC, CT, etc..) into single letter IUPAC codes.
rm(list=ls())
library(seqinr) # bma() Computing an IUPAC nucleotide symbol
inds <- readLines("data/ind393_clean")
dat <- read.delim("data/zos.393ind.HWE.99.geno.code.beagle.gz") # the three genotypes
#dat <- dat[,-(1:4)]

#### make calls
source('R/beagle_calls.fromAllanStrand.R') ### allan's caller. 
beagle_gl(fn="data/zos.393ind.HWE.99.gl.beagle.gz",bamlist= "data/ind393_clean",rowstart=0,nrows=-1,support=log(5)) -> Acalls
loci <- Acalls[,1]
Acalls.tr <- data.frame(t(Acalls[,-1]))

calls <- read.delim('data/zos.393ind.HWE.99.gl.beagle.gz',header=T)
#calls <- calls[,-dim(calls)[2]] # last column is NA
calls[calls=="NN"] <- NA
nind <- 393


for (i in )
out <- c()
  for (j in 1:dim(dat)[1])
  {
    out <- c(out,bma(s2c(dat[j,i])))
  }
write.fasta(out,names = inds[i],file=paste("data/fasta/",inds[i],".fasta",sep=""))




