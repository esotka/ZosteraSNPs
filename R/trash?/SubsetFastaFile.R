######################################
##### subset fasta file of adults #####
######################################
library(ape)
rm(list=ls())
f <- read.FASTA("data/fasta/393ind.fas")
ids <- names(f)
write.FASTA(f[substr(ids,7,7)=="A"],"data/fasta/adults_294ind.fas")
# then to the samples with > 40% basepair calls
### count number of Ns
out <- c()
seqs <- readLines('data/fasta/adults_294ind.fas')[seq(2,2*294,2)]
for (i in 1:294) # 2nd row from fasta file
{
out <- c(out, length(gregexpr(pattern = "N",
                              text=seqs[i])[[1]]))
}

text=read.table('data/fasta/adults_294ind.fas',nrows = 1,skip=1)
n <- nchar(text)

num.bp <- (n-out)/n # the number of basepair calls (i.e., not N)

library(lattice)
histogram(~num.bp,breaks = 50,col="grey")

# list of individuals

ind <- readLines('data/fasta/adults_294ind.fas')[seq(1,2*294-1,2)]

#3 those with less than 10%

alldat <- read.FASTA('data/fasta/adults_294ind.fas')
dat.subset <- alldat[(num.bp)>.4]
write.FASTA(dat.subset,"data/fasta/adults_242ind_foriQtree.fas")
######################################
##### subset fasta file of seeds #####
######################################
rm(list=ls())
f <- read.FASTA("data/fasta/393ind.fas")
ids <- names(f)
write.FASTA(f[substr(ids,7,7)=="S"],"data/fasta/seeds_99ind.fas")
# then to the samples with > 40% basepair calls
### count number of Ns
out <- c()
seqs <- readLines('data/fasta/seeds_99ind.fas')[seq(2,2*99,2)]
for (i in 1:99) # 2nd row from fasta file
{
  out <- c(out, length(gregexpr(pattern = "N",
                                text=seqs[i])[[1]]))
}

text=read.table('data/fasta/seeds_99ind.fas',nrows = 1,skip=1)
n <- nchar(text)

num.bp <- (n-out)/n # the number of basepair calls (i.e., not N)

histogram(~num.bp,breaks = 50,col="grey")

# list of individuals

ind <- readLines('data/fasta/seeds_99ind.fas')[seq(1,2*99-1,2)]

#3 those with less than 10%

alldat <- read.FASTA('data/fasta/seeds_99ind.fas')
dat.subset <- alldat[(num.bp)>.4]
write.FASTA(dat.subset,"data/fasta/seeds_95ind_foriQtree.fas")

