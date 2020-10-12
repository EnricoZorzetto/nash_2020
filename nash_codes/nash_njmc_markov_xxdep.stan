  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years
  int<lower = 1> Nt; // number of years
  real<lower=0> y[M, Nt];  // all events (nyears*nobs)
  // real<lower=-1, upper=1> CORR[M];  // correlations
  real XN[M];  // number of events - sample
  real YN[M];  // number of events - sample
  real XNgen[Mgen];  // number of events - sample
  real YNgen[Mgen];  // number of events - sample
}
parameters {
  // real<lower=0, upper=1> p11;
  // real<lower=0, upper=1> p01;
  // real<lower=0, upper=1> pwet;
  // real<lower=-1, upper=1> corr;
  real pw0;
  real pwx;
  real cr0;
  real crx;
}
transformed parameters{
  // real<lower=0, upper=1> pwet;
  // pwet = p01/(1 - p11 + p01);
  real<lower=0, upper=1> p11[M];
  real<lower=0, upper=1> p01[M];
  real<lower=0, upper=1> pwet[M];
  // real<lower=-1, upper=1> corr[M];
  real<lower=0, upper=1> corr[M];
  real etapw[M];
  real etacr[M];
  
  for (m in 1:M){
    etapw[m] = pw0 + pwx*XN[m];
    pwet[m] =  1/(1+exp(-etapw[m])); // transform between 0 and 1
    
    etacr[m] = cr0 + crx*XN[m];
    // corr[m] =  2/(1+exp(-etacr[m]))-1; // transform between -1 and 1
    corr[m] =  1/(1+exp(-etacr[m])); // transform between 0 and 1
    
    p01[m] = pwet[m]*(1-corr[m]);
    p11[m] = corr[m] + pwet[m]*(1-corr[m]);

 }
}
model {
  // Priors
  // cr0 ~ beta(4, 4);
  // pr ~ beta(4, 4);
  pw0 ~ normal(0, 5);
  pwx ~ normal(0, 1);
  cr0 ~ normal(0, 1);
  crx ~ normal(0, 1);
  
  for (m in 1:M) {
      // first observation (skip for now?)
    if (y[m,1] < 1E-6){
      target += bernoulli_lpmf(0 | pwet[m]);  
    } else {
      target += bernoulli_lpmf(1 | pwet[m]);  
    }
    
    for (j in 2:Nt) {
      if (y[m, j-1] < 1E-6){ // previous is a zero
          if (y[m,j] < 1E-6){ // now ia a zero
            target += bernoulli_lpmf(0 | p01[m]);  
          } else {
            target += bernoulli_lpmf(1 | p01[m]);  
          }
      } else { // previous is a one
          if (y[m,j] < 1E-6){ // mow is a zero
            target += bernoulli_lpmf(0 | p11[m]);  
          } else {
            target += bernoulli_lpmf(1 | p11[m]);  
          }
      }
    }
  }
}
generated quantities {
  vector[M] log_lik;
  // real<lower=-1, upper=1> corr;
  
  
  // corr = p11 - p01;
  
  // vector[Mgen] WETFrep;
  // vector[Mgen] CORRrep;
  // for (m in 1:Mgen) {
  //   WETFrep[m] = bernoulli_rng(pr);
  //   CORRrep[m] = 2*(bernoulli_rng(cr0)-1);
  // }  
  // for (m in 1:M) {
  //   log_lik[m] = bernoulli_lpmf(WETF[m] | pr);
  //   log_lik[m] = log_lik[m] + bernoulli_lpmf(CORR[m] | cr0);
  // }
  
  

    for (m in 1:M) {
    if (y[m,1] < 1E-6){
      log_lik[m] = bernoulli_lpmf(0 | pwet[m]);
    } else {
      log_lik[m] = bernoulli_lpmf(1 | pwet[m]);
    }

    for (j in 2:Nt) {
      if (y[m, j-1] < 1E-6){ // previous is a zero
          if (y[m,j] < 1E-6){ // now ia a zero
            log_lik[m] += bernoulli_lpmf(0 | p01[m]);
          } else {
            log_lik[m] += bernoulli_lpmf(1 | p01[m]);
          }
      } else { // previous is a one
          if (y[m,j] < 1E-6){ // mow is a zero
            log_lik[m] += bernoulli_lpmf(0 | p11[m]);
          } else {
            log_lik[m] += bernoulli_lpmf(1 | p11[m]);
          }
      }
    }
  }

}





