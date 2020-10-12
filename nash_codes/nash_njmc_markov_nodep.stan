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
  real<lower=0, upper=1> pwet;
  // real<lower=0, upper=1> corr01;
  real<lower=0, upper=1> corr;
  // real<lower=-1, upper=1> corr;
}
transformed parameters{
  // real<lower=0, upper=1> pwet;
  // pwet = p01/(1 - p11 + p01);
  // real<lower=-1, upper=1> corr;
  real<lower=0, upper=1> p11;
  real<lower=0, upper=1> p01;
  // corr = 2*corr01-1;
  p01 = pwet*(1-corr);
  p11 = corr + pwet*(1-corr);
  // corr01 = 1+0.5*corr;  // tranformed to be between 0 and 1
}
model {
  // Priors
  pwet ~ beta(4, 4);
  corr ~ beta(2, 10);
  // corr01 ~ beta(20, 20);
  // corr ~ beta(20, 20);
  for (m in 1:M) {
      // first observation (skip for now?)
    if (y[m,1] < 1E-6){
      target += bernoulli_lpmf(0 | pwet);  
    } else {
      target += bernoulli_lpmf(1 | pwet);  
    }
    
    for (j in 2:Nt) {
      if (y[m, j-1] < 1E-6){ // previous is a zero
          if (y[m,j] < 1E-6){ // now ia a zero
            target += bernoulli_lpmf(0 | p01);  
          } else {
            target += bernoulli_lpmf(1 | p01);  
          }
      } else { // previous is a one
          if (y[m,j] < 1E-6){ // mow is a zero
            target += bernoulli_lpmf(0 | p11);  
          } else {
            target += bernoulli_lpmf(1 | p11);  
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
      log_lik[m] = bernoulli_lpmf(0 | pwet);
    } else {
      log_lik[m] = bernoulli_lpmf(1 | pwet);
    }
    
    for (j in 2:Nt) {
      if (y[m, j-1] < 1E-6){ // previous is a zero
          if (y[m,j] < 1E-6){ // now ia a zero
            log_lik[m] += bernoulli_lpmf(0 | p01);
          } else {
            log_lik[m] += bernoulli_lpmf(1 | p01);
          }
      } else { // previous is a one
          if (y[m,j] < 1E-6){ // mow is a zero
            log_lik[m] += bernoulli_lpmf(0 | p11);
          } else {
            log_lik[m] += bernoulli_lpmf(1 | p11);
          }
      }
    }
  }
  
}





