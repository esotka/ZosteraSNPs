### calculation of sigma (σ) = standard deviation of parental position relative to offspring position
# otherwise known as the standard deviation of the dispersal kernel (following Pinsky et al 2016 and Siegel et al 2003. 
rm(list=ls())
niter = 100000
m <- 0.0025465 # isolation by distance slope (Fst/(1-Fst) per km)
mse <- 0.0003331 # standard error on m
a = 0.04 # length of each pop from which seeds were taken 
n = 47 # number of samples generating Ne Estimate - West seeds
NeEstimates = read.table("data/NeEstimator.West.out.txt",header=T,sep="\t") # Nb from NeEstimator with variable MAF

# MAF	Ne	sterr
# 0.05	106.346664	2.254003318
# 0.02	148.1011927	3.377334633
# 0.01	180.8549279	4.613100745
# 0	    180.8549279	4.613100745

# Step 1: De = Effective density per unit km of pop (Ne / km)		
De = NeEstimates$Ne/a
De_se = NeEstimates$sterr/a
# Step 2: estimate sigma
sigma = sqrt(1/(4*De*m))


# two kinds of errors: estimate of Ne (see sterr) and estimate of IBD slope (see mse)
# Following Pinsky et al. 2016, we propagated error through estimations of Ne and m. We used a normal distribution for both.
# used 1 dimensional estimates of sigma
################ FUNCTION 1: generate De  ################
out = c()
for (i in 1:4)
{
De_dist = rnorm(niter, mean = De[i], sd = De_se[i])

################ FUNCTION 2: Generate sigma ################

m_dist = rnorm(niter, mean = m, sd = mse)
sigma_point <- sqrt(1/(4*De[i]*m)) # the point estimate for 1D
sigmas = sqrt(1/(4*De_dist*m_dist)) # the range of estimates for 1D

sigmaU95 = as.numeric(signif(quantile(sigmas,0.975,na.rm=T),3))
sigmaL95 = as.numeric(signif(quantile(sigmas,0.025,na.rm=T),3))


out = rbind(out,data.frame(MAF=NeEstimates$MAF[i],sigma_point, sigmaL95,sigmaU95))
}

print(out)
write.table(out,"output/sigma_estimates.txt",sep="\t",quote=F,row.names=F)

## THESE WERE USING THE NE (ADJUSTED VIA Waples et al. 2014 Genetics)

# n = 214 adults; 6611 SNPs
#   MAF sigma_point sigmaL95 sigmaU95
#1 0.05   0.1921614    0.171    0.223
#2 0.02   0.1628354    0.145    0.189
#3 0.01   0.1473544    0.131    0.171
#4 0.00   0.1473544    0.131    0.171

