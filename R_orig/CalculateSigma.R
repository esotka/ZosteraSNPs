### calculation of sigma (σ) = standard deviation of parental position relative to offspring position
# otherwise known as the standard deviation of the dispersal kernel (following Pinsky et al 2016 and Siegel et al 2003. 
rm(list=ls())
niter = 100000
m <- 0.003151 # isolation by distance slope (Fst/(1-Fst) per km)
mse <- 0.000596 # standard error on m
a = 0.05 # length of each pop from which seeds were taken 
n = 99 # number of samples generating Ne Estimate
NeEstimates = read.table("data/Ne_estimates.txt",header=T,sep="\t") # Nb from NeEstimator with variable MAF

#MAF	Ne	        sterr
#0.05	92.37255068	1.60119027
#0.02	124.7143851	1.445707355
#0.01	161.9074947	1.836600057
#0  	184.6383123	1.98114489


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

### 200 m linear distance
# THESE WERE FROM Nb estimates (unadjusted via Waples et al. 2014 Genetics)
# MAF sigma_point sigmaL95 sigmaU95
# 0.05   0.5117107    0.436    0.645
# 0.02   0.4404374    0.375    0.556
# 0.01   0.3867258    0.330    0.488
# 0.00   0.3621330    0.309    0.456

## THESE WERE USING THE NE (ADJUSTED VIA Waples et al. 2014 Genetics)
# MAF sigma_point sigmaL95 sigmaU95
#1 0.05   0.4144664    0.354    0.524
#2 0.02   0.3566995    0.304    0.449
#3 0.01   0.3130598    0.267    0.394
#4 0.00   0.2931567    0.250    0.369

### 50 m linear distance
#MAF sigma_point sigmaL95 sigmaU95
#1 0.05   0.2072332    0.177    0.261
#2 0.02   0.1783498    0.152    0.225
#3 0.01   0.1565299    0.134    0.197
#4 0.00   0.1465784    0.125    0.185

