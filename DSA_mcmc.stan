// STAN function for fitting data from the manuscript 
//  Analysis of Individual-level Epidemic Data: Study of 2018-2020 Ebola Outbreak in DRC
// by Harley Vossler et al.


// ###################  SIR  ODE Function  

functions {

real[] SIR(real t , real[] y , real[] parms , real[] rdata , int[] idata){

real beta = parms[1];
real gamma = parms[2];
real rho = parms[3];
real dydt[3];
  
 dydt[1] = -beta*y[1]*y[2];
 dydt[2] = beta*y[1]*y[2]-gamma*y[2];
 dydt[3] = gamma*y[2];


return dydt;
 }
}


// ################### Data Module 

data {

  int<lower=0> k;   //number of infected individuals 
  real<lower=0> t0;     // initial time 
  real<lower=0> ti[k];  // times of infection 
  real<lower=0> Ti[k];    // times of recovery 
  int<lower=0,upper=1> event1[k];   // ti censoring variable  (0 if censored 1 ow) 
   int<lower=0,upper=1> event2[k]; //  Ti censoring variable  (0 if censored 1 ow)
}

transformed data{

  real x_r[0];
  int x_i[0];

}

parameters {
  real<lower=0.15> beta;               
  real<lower=0, upper=beta> gamma;               
  real<lower=0, upper=1> rho;   
}

transformed parameters{
  real R0 = beta/gamma;
  real N=1/rho;
  real c = beta*rho;
}

// ##################### Model Defining Module  

model {

  real temp[k,3];
  real parms[3];
  real init[3];
  //real s[k,1];
  //real smax;
  //real factor;


  parms[1] = beta; 
  parms[2]=gamma;
  parms[3]= rho;
  init[1] = 1;
  init[2] = rho; 
  init[3] = 0;

	// ## likelihood function based on SIR ODE 

  temp = integrate_ode_rk45(SIR,init,t0,ti,parms,x_r,x_i,1.0E-7,1.0E-7,1.0E7); //-6

  for(i in 1:k){
  target += event1[i]*(log(temp[i,1])+log(temp[i,2])+log(beta))+event2[i]*log(gamma)-gamma*(Ti[i]-ti[i])+
  (1-event1[i])*(log(temp[i,2]-rho*exp(-gamma*ti[i])));
  
}
 target += -k*log(1-temp[k,1]) ;
 target += gamma_lpdf(beta|.002,.002)+gamma_lpdf(gamma|.002,.002); 

// ## non-informative priors 
}
