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
     col=cols[dat$hab],ylim=c(0,.2))
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

print(summary(m_all))
# y = fst
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.0614417  0.0079542   7.724 3.39e-08 ***
#  x           0.0025876  0.0004857   5.327 1.42e-05 ***
print(summary(m_all_std))
# y = fst/(1-fst)
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept) 0.065455   0.009760   6.706 4.08e-07 ***
#x           0.003151   0.000596   5.286 1.58e-05 ***
print(summary(m_ddss))
# y = fst/(1-fst)
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept) 0.066503   0.019149   3.473  0.00599 **
#  x           0.003089   0.001083   2.854  0.01714 
# y = fst
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept) 0.062839   0.015477   4.060  0.00229 **
#  x           0.002514   0.000875   2.873  0.01657 * 

#######################
##### mantel test #####
#######################
library(mantel)
# all pops
write.csv(dat,"output/Fst~Distance.csv",quote=F,row.names=F)
## convert to matrix manually # this one is for a supplemental table
out = read.table("output/Fst~Distance_matrix.txt",header=T,sep="\t")
rownames(out) = out[,1]
out = out[,-1]

km_mat = as.dist(t(out)) # upper triangle
fst_mat = as.dist(out) # lower triangle
print(mantel(fst_mat,km_mat,permutations=10000))
#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = fst_mat, ydis = km_mat, permutations = 10000) 

#Mantel statistic r: 0.7224 
#      Significance: 0.0012999 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.261 0.375 0.484 0.597 
#Permutation: free
#Number of permutations: 10000

# partial mantel test. Coded habitat at 1 or 0 (different or same habitats)
env = read.table("output/DeepVsShallow_distanceMatrix.txt",header=T,sep="\t")
rownames(env) = env[,1]
env = env[,-1]
env_mat = as.dist(env) # lower triangle
mantel.partial(fst_mat,km_mat,env_mat,permutations = 10000)
mantel.partial(fst_mat,env_mat,km_mat,permutations = 10000)
