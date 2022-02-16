################################################################################
# This script is to accompany the manuscript 
# Analysis of Individual-level Epidemic Data: Study of 2018-2020 Ebola Outbreak in DRC
# by Harley Vossler, Pierre Akilimali, Yuhan Pan, Wasiur R. KhudaBukhsh, 
# Eben Kenah and Grzegorz A. Rempala
################################################################################
#It  fits a single wave of DRC Ebola data via Hamiltonian MCMC using STAN language  
# and saves the MCMS runs.
# 
# Specify the appropriate data paths and source to run.
#
# Written By: Harley Vossler and Grzegorz Rempala 
# Last Edited: Feb 15, 2022
#
################################################################################

### Needed libraries ###

require(rstan); 								# Stan MCMC package
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(1247); 								# Fixed seed for reproducibility

################################################################################

### Input and format data 

file.name= "EVD_Wave3.csv"

dat <- read.csv(file.name);       					
dat = dat[order(dat[,1]),]
T = max(dat[,1],na.rm=T)
k = length(dat[,1])
dat$hosp = pmin(T,dat$hosp)
dat$event2[which(T==dat$hosp)]<-0 				# Check if data needs censoring 

data_SIR <- list(k=k,t0=0,ti=dat$onset,Ti=dat$hosp,event1=dat$event1,event2=dat$event2)

##################################################################################

##############################################################################

### MCMC fit using  STAN library 

fit <- stan(
   file = "DSA_mcmc.stan",  # Stan program
   data = data_SIR,        # named list of data
   chains = 2,             # number of Markov chains
   warmup = 1000,          # number of warmup iterations per chain
   iter = 3000,            # total number of iterations per chain
   cores = 2,              # number of cores
   refresh = 1000,         # show progress every 'refresh' iterations
   control=list(adapt_delta=.9) )

#############################################################################

### output to device  

print(summary(fit))           # output fit summary stats 

                              # output  posterior densities 
print(plot(fit, show_density=T,pars=c("beta","gamma"),
           include=TRUE,fill_color="green")) 

### save final samples to file

out.fit = as.data.frame(cbind(extract(fit)$gamma,
      extract(fit)$beta,extract(fit)$rho,extract(fit)$lp__ ))
colnames(out.fit) = c("gamma","beta","rho","lp")


write.csv(out.fit,paste0("MCMC_Samples_",file.name),row.names=F,col.names=NULL)

#### end script code 