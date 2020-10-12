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
  real c0;
  real cy;
  real cx;
  real w0;
  real wy;
  real wx;
  real n0;
  real ny;
  real nx;
}
transformed parameters{
 real<lower=0> C[M];
 real<lower=0> w[M];
  real<lower=0, upper=1> pn[M];
 real eta[M];

 for(m in 1:M){
 C[m] = exp(c0 + cy*YN[m]*1/(1+exp(-cx*XN[m])));
 w[m] = exp(w0 + wy*YN[m]*1/(1+exp(-wx*XN[m])));
    eta[m] = n0 + ny*YN[m]*1/(1+exp(-nx*XN[m]));
   pn[m] = 1/(1+exp(-eta[m]));
 }
}
model {
  
  n0 ~ normal(0, 1);
  c0 ~ normal(1.5, 0.5);
  w0 ~ normal(-0.4, 0.2);
  
  cx ~ normal(0, 1);
  wx ~ normal(0, 1);
  nx ~ normal(0, 1);
  cy ~ normal(0, 0.25);
  wy ~ normal(0, 0.1);
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
           y[m,j] ~ weibull( w[m], C[m]);
    }
    }
}
}
generated quantities {
  vector[M] log_lik_nj;
  // vector[Mgen] Nrep;
  vector[totnwets] log_lik; // log like for the N model only
  int<lower=0> Nrep[Mgen];
  vector[Mgen] pnrep;
  vector[Mgen] Crep;
  vector[Mgen] wrep;
  vector[Mgen] etarep;
  vector[nrepobs] xijrep;
  vector[Mgen] maxrep;
  int<lower=0> count;
  int<lower=0> startm;
  int<lower=0> endm;
  // int<lower=0> totnwetsrep;

  // totnwetsrep = 0;
  for (m in 1:Mgen) {
    wrep[m] = exp(w0 + wy*YNgen[m]*1/(1+exp(-wx*XNgen[m])));
    Crep[m] = exp(c0 + cy*YNgen[m]*1/(1+exp(-cx*XNgen[m])));
    // Crep[m] = C;
    etarep[m] = n0 + ny*YNgen[m]*1/(1+exp(-nx*XNgen[m]));
    pnrep[m] = 1/(1+exp(-etarep[m]));
    Nrep[m] = binomial_rng(Nt, pnrep[m]);
    // totnwetsrep = totnwetsrep + Nrep[m];
    // generate also daily values - wet at the beginning
    // log_lik[m] = binomial_lpmf(rep_vector(N[m],1) | Nt, pnrep[m]);
    for (j in 1:Nt) {
      if ( j <= Nrep[m]){
        xijrep[ (m-1)*Nt + j ] = weibull_rng(wrep[m], Crep[m]);
      } else {
        xijrep[ (m-1)*Nt + j ] = 0;
      }
    }
    // startm = (m-1)*Nt+1;
    // endm = (m-1)*Nt+Nrep[m];
    // maxrep[m] = max(xijrep[startm:endm]);
    startm = (m-1)*Nt+1;
    endm = (m-1)*Nt+Nrep[m];
    if (Nrep[m] >= 1){
    maxrep[m] = max(xijrep[startm:endm]);
    } else {
    maxrep[m] = 0.0;
    }
  }

  count = 0;
  for (m in 1:M) {
    log_lik_nj[m] = binomial_lpmf(N[m] | Nt, pn[m]);
    for (j in 1:Nt) {
      if (y[m,j] > 1e-6){
        count = count + 1;
        log_lik[count] = weibull_lpdf(y[m,j] | w[m], C[m]);
        }
  }
}
}





