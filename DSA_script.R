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
# Last Edited: Nov 2, 2021
#
################################################################################

### Needed libraries ###

library(rstan)                    # Stan MCMC package
rstan_options(auto_write = TRUE)
set.seed(1247)                    # Fixed seed for repoducibility                         

################################################################################

### Input and format data 

file.name= "Ebola_Wave3.csv"
							
data<-read.csv(file.name);       
data=data[,-1]							
data=data[order(data[,1]),]
T=max(data[,1])
k=length(data[,1])
not.censored=ifelse(data[,2]<=T,1,0) # Check for data censoring  
data[,2]=pmin(T,data[,2])

data_SIR<-list(k=k,t0=0,ti=data[,1],Ti=data[,2],event=not.censored)

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