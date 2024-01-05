### Fst with distance - no Clone adults
rm(list=ls())
library(geosphere) # geographic distance distm()
library(scales) # alpha

# phist (from amova_noClones.R)
dat <- read.csv("output/amova.method1-calls_noClones-pairwisePhiSt.csv")

# get pairwise distance
pos <- data.frame(matrix(c(
  "Curlew_D","-70.91733","42.41907",
  "Curlew_S","-70.91575","42.42009",
  "Lynch_D","-70.85693","42.54366",
  "Lynch_S","-70.85797","42.54488",
  "West_D","-70.80447","42.55789",
  "West_S","-70.80600","42.55922",
  "Niles_D","-70.65666","42.59595",
  "Niles_S","-70.65534","42.59665"
  ),nrow=8,ncol=3,byrow = T))
colnames(pos) <- c("site", "lon","lat")
geo <- distm(pos[,c("lon","lat")])
colnames(geo) <- pos$site; rownames(geo) <- pos$site

dat$km <- c()
for (i in 1:dim(dat)[1])
{dat$km[i] <- geo[match(dat$pop1[i],rownames(geo)),match(dat$pop2[i],colnames(geo))]/1000}



h1 <- unlist(lapply(strsplit(dat$pop1,"_"),"[[",2))
h2 <- unlist(lapply(strsplit(dat$pop2,"_"),"[[",2))
dat$hab <- factor(paste(h1,h2,sep=""))

cols <- c("red","black","pink")
pdf("output/Fst~Distance.pdf",width=4,height=5)
y = dat$phiST; x = dat$km
plot(dat$km,y,xlab="Geographic distance (km)",
     ylab="PHIst/(1-PHIst)",pch=19,col=cols[dat$hab],ylim=c(0,.2))
m = lm(y~x)
#abline(m)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)])
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("lightgrey",0.5),border=alpha("lightgrey",0.5))
m_all = m

par()
y = dat$phiST/(1-dat$phiST); x = dat$km
plot(dat$km,y,xlab="Geographic distance (km)",
     ylab="PHIst/(1-PHIst)",pch=19,col=cols[dat$hab],ylim=c(0,.2))
m = lm(y~x)
#abline(m)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)])
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("lightgrey",0.5),border=alpha("lightgrey",0.5))
m_all_std = m


## shallow only (SS: blue)
ss = dat[dat$hab=="SS",]
y = ss$phiST; x = ss$km
m = lm(y~x)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)],col="blue")
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)],col="blue")
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)],col="blue")
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("blue",0.1),border=alpha("blue",0.1))
m_ss = m

## deep only (DD; red)
dd = dat[dat$hab=="DD",]
y = dd$phiST; x = dd$km
m = lm(y~x)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)],col="red")
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)],col="red")
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)],col="red")
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("red",0.1),border=alpha("red",0.1))
m_dd = m

#### deep-deep and shallow-shallow only (remove deep-vs-shallow)
par()
y = dat$phiST; x = dat$km
plot(dat$km,y,xlab="Geographic distance (km)",
     ylab=expression(paste(Phi," / (1 -",Phi,")")),pch=19,
     col=cols[dat$hab],ylim=c(0,.15))
m = lm(y~x)
#abline(m)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)])
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("lightgrey",0.5),border=alpha("lightgrey",0.5))

ddss = dat[dat$hab%in%c("DD","SS"),]
y = ddss$phiST; x = ddss$km
m = lm(y~x)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)],col="brown")

lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)],col="brown")
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)],col="brown")
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("brown",0.1),border=alpha("brown",0.1))
m_ddss = m



dev.off()
sink("output/Fst~Distance_stats.txt")

print("all data; unstandardized y = fst")
print(summary(m_all))
print("all data; standardized y = fst/(1-fst)")
print(summary(m_all_std))
print("deep-deep and shallow-shallow; unstandardized y = fst")
print(summary(m_ddss))

#######################
##### mantel test #####
#######################
library(ecodist)
# all pops
write.csv(dat,"output/Fst~Distance.csv",quote=F,row.names=F)
## convert to matrix manually # this one is for a supplemental table
out = read.table("output/Fst~Distance_matrix.txt",header=T,sep="\t")
rownames(out) = out[,1]
out = out[,-1]

km_mat = as.dist(t(out)) # upper triangle
fst_mat = as.dist(out) # lower triangle
print("Mantel test - pearson")
print(mantel(fst_mat~km_mat,nperm=10000))
# mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
 # 0.8346216  0.0001000  1.0000000  0.0001000  0.7807977  0.9052253 

# DEFAULT = Pearson correlation

#mantelr	
#Mantel coefficient.

#pval1	
#one-tailed p-value (null hypothesis: r <= 0).

#pval2	
#one-tailed p-value (null hypothesis: r >= 0).

#pval3	
#two-tailed p-value (null hypothesis: r = 0).

#llim	
#lower confidence limit.
#
#ulim	
#upper confidence limit.

sink()