### convert vcf of genotype calls from bcftools view (CC, CT, etc..) into single letter IUPAC codes and then make a fasta file
rm(list=ls())
library(vcfR) # read.vcfR()
library(seqinr) # bma() Computing an IUPAC nucleotide symbol
library(ape)
dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=T)
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

inds <- readLines("data/ind393_clean")
loci.to.use <- readLines("data/loci8684")
#loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")

### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use,]

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

# now pull out individuals with 20% missing data
dna <- read.table("data/fasta/allLoci.txt")

# now make a fasta file
prop <- c()
for (i in 1:dim(dna)[1])
{prop <- c(prop,length(gregexpr("n",dna[i,])[[1]])/nchar(dna[i,]))} # prop of missing data
table(prop<=0.2)
#> table(prop<=0.2)
#FALSE  TRUE 
#66   327

ind20 <- inds[prop<=0.2]
#> table(substr(ind20,7,7))
#A   S 
#232  95 

dna20 <- dna[prop<=0.2,]
write.table(x = dna20,"data/fasta/allLoci_ind20prop.txt",quote=F,row.names = F,col.names = F,sep = "")

for (i in 1:length(ind20))
{
  tmp <- read.table("data/fasta/allLoci_ind20prop.txt",skip=i-1,nrows=1)
  if (i == 1)
  {
    write.table(x = paste(">",ind20[i],"\n",tmp,"\n",sep=""),"data/fasta/327ind8683loc.fas",quote=F,row.names = F,col.names = F)
  }
  else
  {
    write.table(x = paste(">",ind20[i],"\n",tmp,"\n",sep=""),"data/fasta/327ind8683loc.fas",quote=F,row.names = F,col.names = F,append = T)
  }
}

# subset adults and seeds

all <- read.dna("data/fasta/327ind8683loc.fas","fasta")
ad <- all[substr(labels(all),7,7)=="A",]
noClones <- readLines("data/ind245adults_noClones_clean")
ad2 <- ad[labels(ad)%in%noClones,]
write.dna(ad2,"data/fasta/184ind8683loc_noClones_foriQtree.fas","fasta")
seeds <- all[substr(labels(all),7,7)=="S",]
write.dna(seeds,"data/fasta/95seeds_8683loc_foriQtree.fas","fasta")
