###
### take a genotype likelihood file in beagle format and make calls.
### Assumes beagle output created from angsd
### From allan strand - 17 july 2019
library(parallel)
beagle_gl <- function(fn,bamlist=ind,support=log(5),rowstart=0,nrows=-1) #support is the difference in loglike that allows a reliable call.
{
  
  beagle_ind <- function(mal,mil,gl,support)
  {
    names(gl) <- c("h1","het","h2")
    gl=unlist(gl)
    gsup=diff(sort(gl,decreasing=T))[1]
    
    if (abs(gsup)<support)
    {
      NA
    } else {
      
      switch(names(gl)[gl==max(gl)],
             "h1" = paste0(mal,mal),
             "het"= paste0(mal,mil),
             "h2" = paste0(mil,mil))
    }
  }
  
  beagle_snp <- function(r,  angsd_beagle_map = c(A=0,G=2,C=1,T=3),support)
  {
    
    r=unlist(r)
    #        print(r[1:3])
    mal = names(angsd_beagle_map)[r[2]==angsd_beagle_map ]
    mil = names(angsd_beagle_map)[r[3]==angsd_beagle_map ]
    Llikes=log(r[-1:-3])
    starts=seq(1,(length(Llikes)),by=3)
    genos=sapply(starts,function(st){beagle_ind(mal,mil,Llikes[st:(st+2)],support)})
    
    matrix(genos,nrow=1)
  }
  
  rgl =  read.table(gzfile(fn),skip=1+rowstart,nrows=nrows)
  bl = readLines(bamlist)
  inds = as.numeric(gsub("^.*Ap_([0-9]*).*","\\1",bl))
  gmat=as.data.frame(do.call(rbind,mclapply(1:dim(rgl)[1],mc.cores=4,function(i) {beagle_snp(rgl[i,],support=support)})))
  names(gmat) <- paste0("fp",inds)
  retdf=cbind(rgl[,1],gmat)
  names(retdf)[1] ="snp"
  retdf
}
