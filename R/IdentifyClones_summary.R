#### identify clones
rm(list=ls())
### from poppr. see PairwiseDist_identifyClones.r and cloneProcess.R
po = read.csv("data/DissimilarityDist_clones_95perc.csv.gz")
po$adult.seed = paste(substr(po$indA,7,7),substr(po$indB,7,7),sep="")
threshold = 0.05
po = po[po$d<0.05,]
### note: this dataset has duplicates
uniq = combn(unique(po$indA),2)
#a = unique(po$indA); b = unique(po$indB)
po2 = c()
for (i in 1:ncol(uniq))
{
 tmp = po[(po$indA==uniq[1,i] & po$indB==uniq[2,i])| (po$indA==uniq[2,i] & po$indB==uniq[1,i]),][1,]
 if(is.na(tmp[1,1])==FALSE)
       {po2 = rbind(po2,tmp) }
}
po2 = po2[complete.cases(po2),]
#print(dim(po))
#[1] 128   4
#> print(dim(po2))
#[1] 64  4
print(table(po2$adult.seed))
# AA AS 
# 60  4
table(substr(po2$indA,1,3),substr(po2$indB,1,3),po2$adult.seed)

### from ngsRelate. see ngsRelate_processFiles_rab.R
ng = read.table("output/clones_r=80perc.txt",sep="\t",header=T)
print(dim(ng))
# 90  7
print(table(ng$adult.seed))
# AA AS 
# 79 11 

table(substr(ng$indA,1,3),substr(ng$indB,1,3),ng$adult.seed)

# How many samples are clones via both or either method? 

all = rbind(ng[,1:3],po2[c(1,2,4)])
all$method = c(rep("ngsRelate",nrow(ng)),rep("dissimilarity",nrow(po2)))

uniq = combn(unique(c(unique(all$indA),unique(all$indB))),2)
#a = unique(po$indA); b = unique(po$indB)
all2 = c()
for (i in 1:ncol(uniq))
{
 tmp = all[(all$indA==uniq[1,i] & all$indB==uniq[2,i])| (all$indA==uniq[2,i] & all$indB==uniq[1,i]),]
 if (nrow(tmp)==2)
       {
        tmp = tmp[1,]
        tmp[,"method"] = "both"
        }
 if (is.na(tmp[1,1])==FALSE)
       {all2 = rbind(all2,tmp   ) }
}
all2 = all2[complete.cases(all2),]

## switch seed to indA
all2$indA2 = ifelse(all2$adult.seed=="AS" & substr(all2$indA,7,7)=="A",all2$indB,all2$indA) 
all2$indB2 = ifelse(all2$adult.seed=="AS" & substr(all2$indA,7,7)=="A",all2$indA,all2$indB) 
all2$indA = all2$indA2
all2$indB = all2$indB2
all2 = all2[,1:4]
print(table(all2$adult.seed,all2$method))
#     both ngsRelate
#  AA   60        19
#  AS    4         7

all2$siteHabA = substr(all2$indA,1,3)
all2$siteHabB = substr(all2$indB,1,3)
write.csv(all2,"output/Clone_analysis.csv",row.names=F)

### choose clonal replicates to remove

### New list of unique genotypes identified with both methods (dissimilarity and NGS Relate's rab)
# manually #
adults = all2[all2$adult.seed=="AA" & all2$method=="both",]
a.alt <- adults$indA#env$ind[match(rab.sub$a,env$quad.name)]
b.alt <- adults$indB#env$ind[match(rab.sub$b,env$quad.name)]

out <- list()
out[[1]] <- c(a.alt[1],b.alt[1])
for (k in 2:length(a.alt))
{
  a.query <- sapply(out, function(y) a.alt[k] %in% y)
  b.query <- sapply(out, function(y) b.alt[k] %in% y)
    if(sum(a.query) > 0)
    {
      tmp <- unique(c(unlist(out[a.query]),b.alt[k]))
      out[[(1:length(out))[a.query]]] <- tmp
    }
  else
  {
    if(sum(b.query)>0)
      {
      tmp <- unique(c(unlist(out[b.query]),a.alt[k]))
      out[[(1:length(out))[b.query]]] <- tmp
    }
    else
    {
      out[[k]] <- c(a.alt[k],b.alt[k])
    }
  }
}

### get rid of zeros
out <- out[lengths(out)>0]
sink("output/ListOfClones.txt")
print(out)
sink()

### made a list of adults with clones pulled. ***Do not change***

#clones.to.remove <- c()
#clones.to.keep <- c()
#for (m in 1:length(out))
#{
#  tmp <- out[[m]] # keep a randomly chosen clone
#  tmp.to.remove <- tmp[(sample(1:length(tmp),size = length(tmp)-1))]
#  clones.to.remove <- c(clones.to.remove,tmp.to.remove)
#  tmp.to.keep <- tmp[!(tmp%in%tmp.to.remove)]
#  clones.to.keep <- c(clones.to.keep,tmp.to.keep)
#}

#allind = readLines("../datasets/268adults99seeds_clean")
### length(clones.to.remove) = 40
### length(clones.to.keep) = 27
### noClones <- allind[!(allind%in%clones.to.remove)]
### writeLines(noClones,"228adults99seeds_noClones_clean")
#> table(substr(noClones,7,7))
#
#  A   S 
# 214  99
