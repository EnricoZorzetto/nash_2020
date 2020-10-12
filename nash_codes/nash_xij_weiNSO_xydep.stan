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
  real n0;
  real ny;
  real nx;
}
transformed parameters{
 real<lower=0> C[M];
 // real<lower=0> w[M];
  real<lower=0, upper=1> pn[M];
 real eta[M];

 for(m in 1:M){
 C[m] = exp(c0 + cy*YN[m]*1/(1+exp(-cx*XN[m])));
 // w[m] = exp(w0 + wy*YN[m]*1/(1+exp(-wx*XN[m])));
    eta[m] = n0 + ny*YN[m]*1/(1+exp(-nx*XN[m]));
   pn[m] = 1/(1+exp(-eta[m]));
 }
}
model {
  
  w ~ gamma(0.7*10, 10);
  n0 ~ normal(0, 1);
  c0 ~ normal(1.5, 0.5);
  // w0 ~ normal(-0.4, 0.2);
  cx ~ normal(0, 1);
  nx ~ normal(0, 1);
  // wx ~ normal(0, 1);
  cy ~ normal(0, 0.25);
  // wy ~ normal(0, 0.1);
  ny ~ normal(0, 0.5);
  // mc0~inv_gamma(10, 70);
  // sc~inv_gamma(10, 10);
  // mw~inv_gamma(10, 7);
  // sw~inv_gamma(10, 7*0.05); # centered in 0.1
  // Priors
//model for the number of events/year
    for (m in 1:M) {
    // C[m] ~ gumbel(mc, sc); 
    // w[m] ~ gumbel(mw, sw); 
    target += binomial_lpmf(N[m] | Nt, pn[m]);
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
  // int<lower=0> startm;
  // // int<lower=0> endm;
  // int<lower=0> endmN;
  // int<lower=0> endmS;
  // int<lower=0> endmO;
  
    int<lower=0> NrepO[Mgen];
  vector[Mgen] pnrepO;
  vector[Mgen] CrepO;
  vector[Mgen] etarepO;
  // vector[nrepobs] xijrepO;
  real xijrepO[Mgen,Nt];
  vector[Mgen] maxrepO;
  
  int<lower=0> NrepS[Mgen];
  vector[Mgen] pnrepS;
  vector[Mgen] CrepS;
  vector[Mgen] etarepS;
  // vector[nrepobs] xijrepS;
  real xijrepS[Mgen,Nt];
  vector[Mgen] maxrepS;
  
  int<lower=0> NrepN[Mgen];
  vector[Mgen] pnrepN;
  vector[Mgen] CrepN;
  vector[Mgen] etarepN;
  // vector[nrepobs] xijrepN;
  real xijrepN[Mgen,Nt];
  vector[Mgen] maxrepN;
  
  // totnwetsrep = 0;
  for (m in 1:Mgen) {
    CrepO[m] = exp(c0 + cy*YNgen[m]*1/(1+exp(-cx*XNgen[m])));
    CrepS[m] = exp(c0 + cy*YNgenS[m]*1/(1+exp(-cx*XNgenS[m])));
    CrepN[m] = exp(c0 + cy*YNgenN[m]*1/(1+exp(-cx*XNgenN[m])));
    // Crep[m] = C;
    // etarep[m] = n0 + ny*YNgen[m]*1/(1+exp(-nx*XNgen[m]));
    etarepO[m] = n0 + ny*YNgen[m]*1/(1+exp(-nx*XNgen[m]));
    etarepS[m] = n0 + ny*YNgenS[m]*1/(1+exp(-nx*XNgenS[m]));
    etarepN[m] = n0 + ny*YNgenN[m]*1/(1+exp(-nx*XNgenN[m]));
    // pnrep[m] = 1/(1+exp(-etarep[m]));
    pnrepO[m] = 1/(1+exp(-etarepO[m]));
    pnrepS[m] = 1/(1+exp(-etarepS[m]));
    pnrepN[m] = 1/(1+exp(-etarepN[m]));
    // Nrep[m] = binomial_rng(Nt, pnrep[m]);
    NrepO[m] = binomial_rng(Nt, pnrepO[m]);
    NrepS[m] = binomial_rng(Nt, pnrepS[m]);
    NrepN[m] = binomial_rng(Nt, pnrepN[m]);
    
    
    for (j in 1:Nt) {
      if ( j <= NrepO[m]){
        xijrepO[m,j] = weibull_rng(w, CrepO[m]);
      } else {
        xijrepO[m,j] = 0.0;
      }
    }
    for (j in 1:Nt) {
      if ( j <= NrepS[m]){
        xijrepS[m,j] = weibull_rng(w, CrepS[m]);
      } else {
        xijrepS[m,j] = 0.0;
      }
    }
    for (j in 1:Nt) {
      if ( j <= NrepN[m]){
        xijrepN[m,j] = weibull_rng(w, CrepN[m]);
      } else {
        xijrepN[m,j] = 0.0;
      }
    }
    maxrepO[m] = max(xijrepO[m,]);
    maxrepS[m] = max(xijrepS[m,]);
    maxrepN[m] = max(xijrepN[m,]);
  }


    // for (j in 1:Nt) {
    //   // if (NrepS[m] > 0){
    //   if ( j <= NrepO[m]){
    //     xijrepO[ (m-1)*Nt + j ] = weibull_rng(w, CrepO[m]);
    //   } else {
    //     xijrepO[ (m-1)*Nt + j ] = 0.0;
    //   }
    // // }
    // }
    // for (j in 1:Nt) {
    //   // if (NrepS[m] > 0){
    //   if ( j <= NrepS[m]){
    //     xijrepS[ (m-1)*Nt + j ] = weibull_rng(w, CrepS[m]);
    //   } else {
    //     xijrepS[ (m-1)*Nt + j ] = 0.0;
    //   }
    // // }
    // }
    // for (j in 1:Nt) {
    //   // if (NrepN[m] > 0){
    //   if ( j <= NrepN[m]){
    //     xijrepN[ (m-1)*Nt + j ] = weibull_rng(w, CrepN[m]);
    //   } else {
    //     xijrepN[ (m-1)*Nt + j ] = 0.0;
    //   }
    // // }
    // }
    // startm = (m-1)*Nt+1;
    // endmN = (m-1)*Nt+NrepN[m];
    // endmS = (m-1)*Nt+NrepS[m];
    // endmO = (m-1)*Nt+NrepO[m];
    // if (NrepN[m] >= 1){
    // maxrepN[m] = max(xijrepN[startm:endmN]);
    // } else {
    // maxrepN[m] = 0.0;
    // }
    // if (NrepS[m] >= 1){
    // maxrepS[m] = max(xijrepS[startm:endmS]);
    // } else {
    // maxrepS[m] = 0.0;
    // }
    // if (NrepO[m] >= 1){
    // maxrepO[m] = max(xijrepO[startm:endmO]);
    // } else {
    // maxrepO[m] = 0.0;
    // }
    // endmS = (m-1)*Nt+NrepS[m];
    // maxrep[m] = max(xijrep[startm:endm]);
    // maxrepS[m] = max(xijrepS[startm:endmS]);
  // }
  
  count = 0;
  for (m in 1:M) {
    log_lik_nj[m] = binomial_lpmf(N[m] | Nt, pn[m]);
    for (j in 1:Nt) {
      if (y[m,j] > 1e-6){
        count = count + 1;
        log_lik[count] = weibull_lpdf(y[m,j] | w, C[m]);
        }
  }
}
}





