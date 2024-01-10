## analysis of spagedi SGS data
## 10 sites
rm(list=ls())
dat <- readLines("data/allSite.Spagedi.out.txt")
dat2 <- strsplit(dat,"\t")
ncat <- 10

# 2-sided test
tmp = dat[grep("2-sided test",dat)]
sitenames = substr(tmp,regexpr("out.txt",tmp)-6,regexpr("out.txt",tmp)-2)
p = c()
for(i in 3:12)
{p <- c(p,as.numeric(unlist(lapply(strsplit(tmp,"\t"),"[[",i))))}
p <- matrix(p,nrow=ncat,ncol=8,byrow=T)
p <- ifelse(p<=0.05,"*","")
colnames(p) = sitenames

# meters
tmp = dat[grep("Mean distance",dat)]
sitenames = substr(tmp,regexpr("out.txt",tmp)-6,regexpr("out.txt",tmp)-2)
m = c()
for(i in 3:12)
{m <- c(m,as.numeric(unlist(lapply(strsplit(tmp,"\t"),"[[",i))))}
m <- matrix(m,nrow=ncat,ncol=8,byrow=T)
colnames(m) = sitenames

# kinship
tmp = dat[grep("ALL LOCI",dat)]
sitenames = substr(tmp,regexpr("out.txt",tmp)-6,regexpr("out.txt",tmp)-2)
k = c()
for(i in 3:12)
{k <- c(k,as.numeric(unlist(lapply(strsplit(tmp,"\t"),"[[",i))))}
k <- matrix(k,nrow=ncat,ncol=8,byrow=T)
colnames(k) = sitenames

# lower CI
tmp = dat[grep("95%CI-inf",dat)]
sitenames = substr(tmp,regexpr("out.txt",tmp)-6,regexpr("out.txt",tmp)-2)
l95 = c()
for(i in 3:12)
{l95 <- c(l95,as.numeric(unlist(lapply(strsplit(tmp,"\t"),"[[",i))))}
l95 <- matrix(l95,nrow=ncat,ncol=8,byrow=T)
colnames(l95) = sitenames

# upper CI
tmp = dat[grep("95%CI-sup",dat)]
sitenames = substr(tmp,regexpr("out.txt",tmp)-6,regexpr("out.txt",tmp)-2)
u95 = c()
for(i in 3:12)
{u95 <- c(u95,as.numeric(unlist(lapply(strsplit(tmp,"\t"),"[[",i))))}
u95 <- matrix(u95,nrow=ncat,ncol=8,byrow=T)
colnames(u95) = sitenames

# pretty names
pretty_names = c("Curlew Deep",#cur_d" 
                 "Curlew Shallow",#,cur_s" 
                 "Lynch Deep",#lyn_d" 
                 "Lynch Shallow",#lyn_s" 
                 "West Deep",#wes_d"
                 "West Shallow",#wes_s" 
                 "Niles Deep",#"nil_d" 
                 "Niles Shallow")#nil_s"

pdf("output/SpagediSGS_noClones_10categories.pdf",width=10,height=4)
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1))
for(n in c(1:4,7,8,5,6)) # change order of sites from South to North
{
  plot(x=m[,n],y=k[,n],type="b",ylim=c(-.1,.15),xlim=c(0,45),ylab="kinship",xlab="meters",pch=20,cex=1.2)
  mtext(pretty_names[n],line=-2)
  segments(-1,0,50,0,lty="dashed")
  points(x=m[,n],y=u95[,n],col="red",lty="dashed",type="l")
  points(x=m[,n],y=l95[,n],col="red",lty="dashed",type="l")
  text(x=m[,n],y=-0.1,p[,n],cex=1.5)
}
dev.off()


