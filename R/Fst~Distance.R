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

cols <- c("red","black","blue")
pdf("output/Fst~Distance.pdf",width=7,height=5)
y = dat$phiST/(1-dat$phiST); x = dat$km
plot(dat$km,y,xlab="Geographic distance (km)",
     ylab="Genetic Distance (Fst/(1-Fst))",pch=19,col=cols[dat$hab],ylim=c(0,.2))
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
     ylab="Genetic Distance (Fst/(1-Fst))",pch=19,col=cols[dat$hab],ylim=c(0,.2))
m = lm(y~x)
#abline(m)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)])
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("lightgrey",0.5),border=alpha("lightgrey",0.5))
m_all = m


## shallow only (SS: blue)
ss = dat[dat$hab=="SS",]
y = ss$phiST/(1-ss$phiST); x = ss$km
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
y = dd$phiST/(1-dd$phiST); x = dd$km
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
y = dat$phiST/(1-dat$phiST); x = dat$km
plot(dat$km,y,xlab="Geographic distance (km)",
     ylab="Genetic Distance (Fst/(1-Fst))",pch=19,col=cols[dat$hab],ylim=c(0,.2))
m = lm(y~x)
#abline(m)
lines(x[order(x)], predict(m, interval = "confidence")[, "fit"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "lwr"][order(x)])
lines(x[order(x)], predict(m, interval = "confidence")[, "upr"][order(x)])
polygon(x = c(x[order(x)],x[order(x,decreasing = T)]),
        y=c(predict(m, interval = "confidence")[, "lwr"][order(x)],predict(m, interval = "confidence")[, "upr"][order(x,decreasing = T)]),
        col=alpha("lightgrey",0.5),border=alpha("lightgrey",0.5))

ddss = dat[dat$hab%in%c("DD","SS"),]
y = ddss$phiST/(1-ddss$phiST); x = ddss$km
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
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.065455   0.009760   6.706 4.08e-07 ***
#x           0.003151   0.000596   5.286 1.58e-05 ***
## Slope CI = mean ± 2*SE 
## > 0.003151 + 2*0.000596
#[1] 0.004343
#> 0.003151 - 2*0.000596
#[1] 0.001959
# Intercept
#> 0.065455 + 2*0.009760
#[1] 0.084975
#> 0.065455 - 2*0.009760
#[1] 0.045935
print(summary(m_ss))
#Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.039746   0.021534   1.846   0.1387  
#x           0.003934   0.001217   3.233   0.0319 *
print(summary(m_dd))
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.093095   0.027376   3.401   0.0273 *
#  x           0.002254   0.001548   1.456   0.2192  
print(summary(m_ddss))
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept) 0.066503   0.019149   3.473  0.00599 **
#  x           0.003089   0.001083   2.854  0.01714 

y = dat$phiST/(1-dat$phiST); x = dat$km
print(cor.test(y,x,method="spearman"))
y = ss$phiST/(1-ss$phiST); x = ss$km
print(cor.test(y,x,method="spearman"))
y = dd$phiST/(1-dd$phiST); x = dd$km
print(cor.test(y,x,method="spearman"))
y = ddss$phiST/(1-ddss$phiST); x = ddss$km
print(cor.test(y,x,method="spearman"))



dat$hab[dat$hab=="SD"] <- "DS"
dat$hab <- factor(dat$hab)
anova(lm(phiST~hab,dat)) # no significant difference between DD, DS, SS

#Analysis of Variance Table

#Response: phiST
#Df    Sum Sq    Mean Sq F value Pr(>F)
#hab        2 0.0021507 0.0010753  1.0599 0.3616
#Residuals 25 0.0253647 0.0010146
anova(lm(phiST~hab,dat[!dat$hab=="DS",])) # no significant difference between DD or SS

#Analysis of Variance Table

#Response: phiST
#Df    Sum Sq    Mean Sq F value Pr(>F)
#hab        1 0.0013659 0.00136589  1.8154 0.2076
#Residuals 10 0.0075239 0.00075239              

## we don't have much power to detect differences in the slope of these (low sample size)

