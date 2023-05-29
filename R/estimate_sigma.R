######################################################################
## Resampling approach to estimate 95% CI on dispersal kernel spread
######################################################################
## Borrowed from Pinsky et al 2016 Current Biology
rm(list=ls())
# Set parameters
niter <- 1000000 # number of iterations
alpha <- 2 # minimum age of maturity (for Waples Nb to Ne calculations)
AL <- 3 # adult lifespan (for Waples Nb to Ne calculations)
A <- .2 # length of a population where seeds were collected (km)
m <- 0.002514 # isolation by distance slope (genetic distance/km)
mse <- 0.000875 # standard error on m
Waples_Nb_est <- 109.2 #(west beach; MAF > 0.01)
Waples_Nbl95_est <- 83.1
Waples_Nbu95_est <- 152.1
#Waples_Nb_est <- 61.8 #(west beach; MAF > 0.05)
#Waples_Nbl95_est <- 42.5
#Waples_Nbu95_est <- 97


########### FUNCTION ChiSqCIs() ###########

findChiSqCIs <- function(De, Del95, Deu95, verbose = FALSE){
a2 = 105
a1 = 100
while (abs(mean(as.numeric(a1))-a2) > 1){ # inefficient search for a df that fits my CIs
  if(mean(as.numeric(a1)) > a2) a2 = a2 - 1
  if(mean(as.numeric(a1)) < a2) a2 = a2 + 1
  a1 = c(Del95, Deu95)*qchisq(c(0.975, 0.025), df=a2)/De
  if(verbose) print(paste(paste(round(mean(as.numeric(a1)),2), collapse=','), a2))
}
return(a2)
}

########### FUNCTION DeFromNb() ###########

# Calculate sigma (dispersal kernel spread) from number of breeders (Nb)
# Uses Waples et al. 2013 PRSB/Waples et al. 2014 Genetics correction for Nb to Ne
#
# niter: number of iterations to produce
# AL: adult lifespan (years)
# alpha: minimum age of maturity (years)
# Waples_Nb: linkage disequilibrium estimate of the number of breeders (e.g., output from NeEstimator based on a single cohort)
# Waples_Nbl95: lower 95% confidence interval for Waples_Nb
# Waples_Nbu95: upper 95% confidence interval for Waples_Nb
# A: length (1D) or area (2D) of the region to which the Nb estimate applies


DeFromNb <- function(niter, AL, alpha, Waples_Nb, Waples_Nbl95, Waples_Nbu95, A){

  # Fit the Waples et al. model to relate log10(AL/alpha) to the Nb/Ne ratio
  dat = read.csv('data/Pinsky/AllSpecies4_trim.csv') # data from Waples et al. 2013 PRSB
  mod = lm(NbNe ~ logRLAM, data=dat) # fit model
  
  # Get CIs and samples for our Nb/Ne ratio (using AL and alpha)
  NbNe = predict(mod, new=data.frame(logRLAM = log10(AL/alpha)), interval='prediction', level=0.95) # get the 95%CI prediction intervals
  a = c(NbNe[,1] - NbNe[,2], NbNe[,3] - NbNe[,1]) # the CI intervals
  NbNes1 = rnorm(niter, mean = NbNe[,1], sd = mean(a)/1.96) # generate first set of Nb/Ne samples
  NbNes2 = rnorm(niter, mean = NbNe[,1], sd = mean(a)/1.96) # generate second set (need two draws since used 2x in equation)
  
  # generate point values for De
  Waples_Nb_adj = Waples_Nb/(1.26-0.323*NbNe[,1]) # point estimate for Nb after correction. From Waples et al. 2014 Genetics. We use the equation for "True Nb/Ne" from Table 3 because we have estimated Nb/Ne as a probability distribution (the "Using two traits" version in Table 3 is nothing more than inserting Nb/Ne=0.485+0.758log(AL/alpha) into Eq. 8). 
  Waples_Ne_adj = Waples_Nb_adj/NbNe[,1] # point estimate for Ne after correction. From Eq. Waples et al. 2014. We used the equation for "True Nb/Ne" from Table 3 because we have estimated Nb/Ne as a probability distribution.
  De <- Waples_Ne_adj/A # the point estimate for De
  
  # generate samples from Nb
  a2 <- findChiSqCIs(Waples_Nb, Waples_Nbl95, Waples_Nbu95, verbose=FALSE) # find the degrees of freedom for a chi-squared distribution that fits our upper and lower 95% CI values
  WaplesNbs = a2*Waples_Nb/rchisq(niter, df=a2) # generate Nb values from chisq distribution (unadjusted Nb)
  
  # optional debugging to check fitting of the chi-squared distribution
  # print(paste('a2=', a2, 'for i=',i))
  # print(rbind(signif(quantile(WaplesNbs, c(0.025, 0.975)),4), signif(Waples_Nbl95, Waples_Nbu95))) # second line values should be close to values in first line
  
  # convert range of Nbs to range of Nes and Des using Waples et al. 2013 PRSB equations
  Nbs_adj = WaplesNbs/(1.26+0.323*NbNes1)
  Nes_adj = Nbs_adj/NbNes2
  Des <- Nes_adj/A # and for De
  #	Del95 = quantile(Des, 0.025)
  #	Deu95 = quantile(Des, 0.975)
  
  return(list(De=De, Des=Des))
  
}

########### FUNCTION sigmaFrom_m() ###########
# Estimate sigma from effective density and isolation by distance slope
#
# De: point estimate of effective density
# Des: samples from the probability distribution of De values
# m: the estimate of the isolation by distance slope
# mse: the standard error on m
# dims: the number of dimensions (1 or 2)


sigmaFrom_m <- function(De, Des=rep(De,1000), m, mse, dims=1){
  niter <- length(Des)
  
  # generate range of slope estimates
  ms = rnorm(niter, mean = m, sd = mse)
  
  # calculate sigma
  if(dims==1){
    sigma_point <- sqrt(1/(4*De*m)) # the point estimate for 1D
    sigmas = sqrt(1/(4*Des*ms)) # the range of estimates for 1D
  }
  if(dims==2){
    sigma_point = sqrt(1/(4*pi*De*m)) # the point estimate for 2D
    sigmas = sqrt(1/(4*pi*Des*ms)) # the range of estimates for 2D
  }
  if(!(dims %in% c(1,2))){
    stop('Dimensions need to be 1 or 2')
  }
  
  return(list(sigma_point=sigma_point, sigmas=sigmas))
}

########### FUNCTION summarizeSigmas() ###########

# Print some useful statistics for a distribution of sigma values
#
# sigmas: the samples from the sigma distribution
# sg: the number of significant digits to report

summarizeSigmas <- function(sigmas, sg=3){
  
  sigma_mean = signif(mean(sigmas, na.rm=TRUE), sg)
  sigma_median = signif(median(sigmas, na.rm=TRUE), sg)
  sigma_sd = signif(sd(sigmas, na.rm=TRUE),3)
  sigma_l95 = signif(quantile(sigmas, 0.025, na.rm=TRUE), sg)
  sigma_u95 = signif(quantile(sigmas, 0.975, na.rm=TRUE), sg)
  
  cat(paste(
    'Mean:         ', sigma_mean, '\n',
    'Median:       ', sigma_median, '\n',
    'SD:           ', sigma_sd, '\n',
    'Lower 95% CI: ', sigma_l95, '\n',
    'Upper 95% CI: ', sigma_u95, sep=''))
  
}




# Estimate sigma (dispersal kernel spread) in various ways

# Using Waples Nb (# breeders) and IBD slope
# Waples_Nb value and CI from NeEstimator results
D1 <- DeFromNb(niter, AL, alpha, Waples_Nb_est, Waples_Nbl95_est, Waples_Nbu95_est, A) # estimate De
s1 <- sigmaFrom_m(De=D1$De, Des=D1$Des, m=m, mse=mse, dims=1) # estimate sigma
print(s1$sigma_point) # report point value for sigma
print(summarizeSigmas(s1$sigmas,3)) # report distribution for sigma
