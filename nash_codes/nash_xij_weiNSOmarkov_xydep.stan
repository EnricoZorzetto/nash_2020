  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of years
  real<lower=0> y[M,Nt]; // all rainfall values
  int<lower=0> N[M];  // number of events - sample
  real XN[M];  
  real YN[M];  
  real XNgen[Mgen];  
  real YNgen[Mgen];  
  real XNgenS[Mgen];  
  real YNgenS[Mgen];  
  real XNgenN[Mgen];  
  real YNgenN[Mgen];  
  int<lower = 1> totnwets; // total number of events above zero in the  
    // prior parameters for hyper-distributions (shape and rate)
    // real<lower = 0> mc0prior[2]; // a, b of an inverse gamma
    // real<lower = 0> sc0prior[2]; // a, b of an inverse gamma
    // real<lower = 0> mw0prior[2]; // a, b of an inverse gamma
    // real<lower = 0> sw0prior[2]; // a, b of an inverse gamma
    // real<lower = 0> pn0prior[2];
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
transformed data{
  int<lower=0> nrepobs;
  nrepobs = Nt*Mgen;
}
parameters {
  real<lower=0> w;
  real c0;
  real cy;
  real cx;
  // real w0;
  // real wy;
  // real wx;
  // real n0;
  // real ny;
  // real nx;
  real pw0;
  real pwy;
  real pwx;
  real cr0;
  real cry;
  real crx;
}
transformed parameters{
 real<lower=0> C[M];
 // real<lower=0> w[M];
  // real<lower=0, upper=1> pn[M];
 // real eta[M];
 
   real<lower=0, upper=1> p11[M];
  real<lower=0, upper=1> p01[M];
  real<lower=0, upper=1> pwet[M];
  // real<lower=-1, upper=1> corr[M];
  real<lower=0, upper=1> corr[M];
  real etapw[M];
  real etacr[M];

 for(m in 1:M){
   
    // model for rainfall intensities
    C[m] = exp(c0 + cy*YN[m]*1/(1+exp(-cx*XN[m])));
    // w[m] = exp(w0 + wy*YN[m]*1/(1+exp(-wx*XN[m])));
    // eta[m] = n0 + ny*YN[m]*1/(1+exp(-nx*XN[m]));
    // pn[m] = 1/(1+exp(-eta[m]));
   
       
    // Markov Chain for rainfall arrivals
    etapw[m] = pw0 + pwy*YN[m]*1/(1+exp(-pwx*XN[m]));
    pwet[m] =  1/(1+exp(-etapw[m])); // transform between 0 and 1
    etacr[m] = cr0 + cry*YN[m]*1/(1+exp(-crx*XN[m]));
    // corr[m] =  2/(1+exp(-etacr[m]))-1; // transform between -1 and 1
    corr[m] =  1/(1+exp(-etacr[m])); // transform between 0 and 1
    p01[m] = pwet[m]*(1-corr[m]);
    p11[m] = corr[m] + pwet[m]*(1-corr[m]);
 }
}
model {
  
  pw0 ~ normal(0, 5);
  pwy ~ normal(0, 1);
  pwx ~ normal(0, 1);
  cr0 ~ normal(0, 1);
  cry ~ normal(0, 1);
  crx ~ normal(0, 1);
  
  w ~ gamma(0.7*10, 10);
  // n0 ~ normal(0, 1);
  c0 ~ normal(1.5, 0.5);
  // w0 ~ normal(-0.4, 0.2);
  cx ~ normal(0, 1);
  // nx ~ normal(0, 1);
  // wx ~ normal(0, 1);
  cy ~ normal(0, 0.25);
  // wy ~ normal(0, 0.1);
  // ny ~ normal(0, 0.5);
  // mc0~inv_gamma(10, 70);
  // sc~inv_gamma(10, 10);
  // mw~inv_gamma(10, 7);
  // sw~inv_gamma(10, 7*0.05); # centered in 0.1
  // Priors
  
//model for the number of events/year
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
  
//model for the intensities
    for (m in 1:M) {
      for (j in 1:Nt) {
        if (y[m,j] >1e-6) {
           y[m,j] ~ weibull( w, C[m]);
    }
    }
}
}
generated quantities {
  vector[M] log_lik_nj;
  // vector[Mgen] Nrep;
  // int<lower=0> Nrep[Mgen];
  vector[totnwets] log_lik; // log like for the N model only
  // vector[Mgen] pnrep;
  // vector[Mgen] Crep;
  // vector[Mgen] wrep;
  // vector[Mgen] etarep;
  // vector[nrepobs] xijrep;
  // vector[Mgen] maxrep;
  int<lower=0> count;
  int<lower=0> token;
  // int<lower=0> startm;
  // // int<lower=0> endm;
  // int<lower=0> endmN;
  // int<lower=0> endmS;
  // int<lower=0> endmO;
  
    int<lower=0> NrepO[Mgen];
  // vector[Mgen] pnrepO;
  vector[Mgen] CrepO;
  // vector[Mgen] etarepO;
  // vector[nrepobs] xijrepO;
  real xijrepO[Mgen, Nt];
  vector[Mgen] maxrepO;
  
  int<lower=0> NrepS[Mgen];
  // vector[Mgen] pnrepS;
  vector[Mgen] CrepS;
  // vector[Mgen] etarepS;
  // vector[nrepobs] xijrepS;
  real xijrepS[Mgen, Nt];
  vector[Mgen] maxrepS;
  
  int<lower=0> NrepN[Mgen];
  // vector[Mgen] pnrepN;
  vector[Mgen] CrepN;
  // vector[Mgen] etarepN;
  // vector[nrepobs] xijrepN;
  real xijrepN[Mgen, Nt];
  vector[Mgen] maxrepN;
  
  real<lower=0, upper=1> p11repO[Mgen];
  real<lower=0, upper=1> p01repO[Mgen];
  real<lower=0, upper=1> pwetrepO[Mgen];
  real<lower=0, upper=1> corrrepO[Mgen];
  real etapwrepO[Mgen];
  real etacrrepO[Mgen];
  
    
  real<lower=0, upper=1> p11repN[Mgen];
  real<lower=0, upper=1> p01repN[Mgen];
  real<lower=0, upper=1> pwetrepN[Mgen];
  real<lower=0, upper=1> corrrepN[Mgen];
  real etapwrepN[Mgen];
  real etacrrepN[Mgen];
  
    
  real<lower=0, upper=1> p11repS[Mgen];
  real<lower=0, upper=1> p01repS[Mgen];
  real<lower=0, upper=1> pwetrepS[Mgen];
  real<lower=0, upper=1> corrrepS[Mgen];
  real etapwrepS[Mgen];
  real etacrrepS[Mgen];
  
  
  

  for (m in 1:Mgen) {
    // generate intensity parameters - ONS
    CrepO[m] = exp(c0 + cy*YNgen[m]*1/(1+exp(-cx*XNgen[m])));
    CrepS[m] = exp(c0 + cy*YNgenS[m]*1/(1+exp(-cx*XNgenS[m])));
    CrepN[m] = exp(c0 + cy*YNgenN[m]*1/(1+exp(-cx*XNgenN[m])));
    
    // generate arrivals parameters - O
    etapwrepO[m] = pw0 + pwy*YNgen[m]*1/(1+exp(-pwx*XNgen[m]));
    pwetrepO[m] =  1/(1+exp(-etapw[m])); // transform between 0 and 1
    etacrrepO[m] = cr0 + cry*YNgen[m]*1/(1+exp(-crx*XNgen[m]));
    corrrepO[m] =  1/(1+exp(-etacrrepO[m])); // transform between 0 and 1
    p01repO[m] = pwetrepO[m]*(1-corrrepO[m]);
    p11repO[m] = corrrepO[m] + pwetrepO[m]*(1-corrrepO[m]);
    
    // generate arrivals parameters - N
    etapwrepS[m] = pw0 + pwy*YNgenS[m]*1/(1+exp(-pwx*XNgenS[m]));
    pwetrepS[m] =  1/(1+exp(-etapw[m])); // transform between 0 and 1
    etacrrepS[m] = cr0 + cry*YNgenS[m]*1/(1+exp(-crx*XNgenS[m]));
    corrrepS[m] =  1/(1+exp(-etacrrepS[m])); // transform between 0 and 1
    p01repS[m] = pwetrepS[m]*(1-corrrepS[m]);
    p11repS[m] = corrrepS[m] + pwetrepS[m]*(1-corrrepS[m]);
    
    // generate arrivals parameters - S
    etapwrepN[m] = pw0 + pwy*YNgenN[m]*1/(1+exp(-pwx*XNgenN[m]));
    pwetrepN[m] =  1/(1+exp(-etapw[m])); // transform between 0 and 1
    etacrrepN[m] = cr0 + cry*YNgenN[m]*1/(1+exp(-crx*XNgenN[m]));
    corrrepN[m] =  1/(1+exp(-etacrrepN[m])); // transform between 0 and 1
    p01repN[m] = pwetrepN[m]*(1-corrrepN[m]);
    p11repN[m] = corrrepN[m] + pwetrepN[m]*(1-corrrepN[m]);
    
    
    
    // // generate first day of the season (wet / dry)
    // NrepO[m] = 0;
    // token = bernoulli_rng(pwet[m]);
    // if (token==1){
    //   NrepO[m] += 1;
    //   xijrepO[m, 1] = weibull_rng(w, CrepO[m]);
    // } else {
    //   xijrepO[m, 1] = 0.0;
    // }
    // for (j in 2:Nt){
    //   if (token==1){ //previous day was wet
    //     token = bernoulli_rng(p11repO[m]); // if 1 it stays 1
    //       if (token==1){
    //         NrepO[m] += 1;
    //         xijrepO[m, j] = weibull_rng(w, CrepO[m]);
    //       } else {
    //         xijrepO[m, j] = 0.0;
    //       }
    //     } else { // previous day dry
    //       token = bernoulli_rng(p01repO[m]); // if O becomes 1
    //       if (token==1){
    //         NrepO[m] += 1;
    //         xijrepO[m, j] = weibull_rng(w, CrepO[m]);
    //       } else {
    //         xijrepO[m, j] = 0.0;
    //       }
    //     }
    //   }
    // maxrepO[m] = max(xijrepO[m,]);  // compute the maxima for this season:
    
    
    
    // generate first day of the season (wet / dry)
    NrepO[m] = 0;
    token = bernoulli_rng(pwetrepO[m]);
    if (token==1){
      NrepO[m] += 1;
      xijrepO[m, 1] = weibull_rng(w, CrepO[m]);
    } else {
      xijrepO[m, 1] = 0.0;
    }
    for (j in 2:Nt){
      if (token==1){ //previous day was wet
        token = bernoulli_rng(p11repO[m]); // if 1 it stays 1
          if (token==1){
            NrepO[m] += 1;
            xijrepO[m, j] = weibull_rng(w, CrepO[m]);
          } else {
            xijrepO[m, j] = 0.0;
          }
        } else { // previous day dry
          token = bernoulli_rng(p01repO[m]); // if O becomes 1
          if (token==1){
            NrepO[m] += 1;
            xijrepO[m, j] = weibull_rng(w, CrepO[m]);
          } else {
            xijrepO[m, j] = 0.0;
          }
        }
      }
    maxrepO[m] = max(xijrepO[m,]);  // compute the maxima for this season:
    
    
        // generate first day of the season (wet / dry)
    NrepS[m] = 0;
    token = bernoulli_rng(pwetrepS[m]);
    if (token==1){
      NrepS[m] += 1;
      xijrepS[m, 1] = weibull_rng(w, CrepS[m]);
    } else {
      xijrepS[m, 1] = 0.0;
    }
    for (j in 2:Nt){
      if (token==1){ //previous day was wet
        token = bernoulli_rng(p11repS[m]); // if 1 it stays 1
          if (token==1){
            NrepS[m] += 1;
            xijrepS[m, j] = weibull_rng(w, CrepS[m]);
          } else {
            xijrepS[m, j] = 0.0;
          }
        } else { // previous day dry
          token = bernoulli_rng(p01repS[m]); // if O becomes 1
          if (token==1){
            NrepS[m] += 1;
            xijrepS[m, j] = weibull_rng(w, CrepS[m]);
          } else {
            xijrepS[m, j] = 0.0;
          }
        }
      }
    maxrepS[m] = max(xijrepS[m,]);  // compute the maxima for this season:
    
    
    
    // generate first day of the season (wet / dry)
    NrepN[m] = 0;
    token = bernoulli_rng(pwetrepN[m]);
    if (token==1){
      NrepN[m] += 1;
      xijrepN[m, 1] = weibull_rng(w, CrepN[m]);
    } else {
      xijrepN[m, 1] = 0.0;
    }
    for (j in 2:Nt){
      if (token==1){ //previous day was wet
        token = bernoulli_rng(p11repN[m]); // if 1 it stays 1
          if (token==1){
            NrepN[m] += 1;
            xijrepN[m, j] = weibull_rng(w, CrepN[m]);
          } else {
            xijrepN[m, j] = 0.0;
          }
        } else { // previous day dry
          token = bernoulli_rng(p01repN[m]); // if O becomes 1
          if (token==1){
            NrepN[m] += 1;
            xijrepN[m, j] = weibull_rng(w, CrepN[m]);
          } else {
            xijrepN[m, j] = 0.0;
          }
        }
      }
    maxrepN[m] = max(xijrepN[m,]);  // compute the maxima for this season:
  
  

    
    
    
  } // end loop on years



  
    // LIKELIHOOD FOR MODEL FOR RAINFALL INTENSITIES
  count = 0;
  for (m in 1:M) {
    for (j in 1:Nt) {
      if (y[m,j] > 1e-6){
        count = count + 1;
        log_lik[count] = weibull_lpdf(y[m,j] | w, C[m]);
        }
  }
}

    // LIKELIHOOD FOR MARKOV CHAIN MODEL
    for (m in 1:M) {
    if (y[m,1] < 1E-6){
      log_lik_nj[m] = bernoulli_lpmf(0 | pwet[m]);
    } else {
      log_lik_nj[m] = bernoulli_lpmf(1 | pwet[m]);
    }

    for (j in 2:Nt) {
      if (y[m, j-1] < 1E-6){ // previous is a zero
          if (y[m,j] < 1E-6){ // now ia a zero
            log_lik_nj[m] = log_lik_nj[m] + bernoulli_lpmf(0 | p01[m]);
          } else {
            log_lik_nj[m] = log_lik_nj[m] + bernoulli_lpmf(1 | p01[m]);
          }
      } else { // previous is a one
          if (y[m,j] < 1E-6){ // mow is a zero
            log_lik_nj[m] = log_lik_nj[m] + bernoulli_lpmf(0 | p11[m]);
          } else {
            log_lik_nj[m] = log_lik_nj[m] + bernoulli_lpmf(1 | p11[m]);
          }
      }
    }
  }
}





