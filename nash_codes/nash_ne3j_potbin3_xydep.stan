  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years
  int<lower = 1> Nt; // number of events /year
  int<lower=0> NP3[M];  // number of events - sample
  real XN[M];  // number of events - sample
  real YN[M];  // number of events - sample
  real XNgen[Mgen];  // number of events - sample
  real YNgen[Mgen];  // number of events - sample
}
parameters {
  real n3x;
  // real sy;
  // real<lower=0, upper=1> mn0;
  real n30;
  real n3y;
  // real<lower=0, upper=1> ma2;
  // real<lower=0> bn;
  // real<lower=0, upper = 1> pn[M];
}
transformed parameters{
 // real<lower=0> an[M];
 // real<lower=0> bn[M];
 real<lower=0, upper=1> pn3[M];
 // real ma;
 real eta3[M];
 // real mnv[M];
 // ma = 0.5*ma2;
 for (m in 1:M){

   // eta[m] = eta0 + etay*YN[m];
   eta3[m] = n30 + n3y*YN[m]*1/(1+exp(-n3x*XN[m]));
   pn3[m] = 1/(1+exp(-eta3[m]));
 }
}
model {
  // Priors
  // mn0 ~ beta(2, 2);
  // ma2 ~ beta(2, 2);
  // bn ~ gamma(10,1); // I want something fairly uninformative
  n3x ~ normal(0, 2);
  // sy ~ normal(0, 1);
  for (m in 1:M) {
    // pn[m] ~ beta(an[m], bn);  
  target += binomial_lpmf(NP3[m] | Nt, pn3[m]);
  }
}
generated quantities {
  vector[M] log_lik;
  vector[Mgen] Nrep;
  vector[Mgen] pnrep;
  vector[Mgen] etarep;
  // vector[Nt] predpmf;
  for (m in 1:Mgen) {
    etarep[m] = n30 + n3y*YNgen[m]*1/(1+exp(-n3x*XNgen[m]));
    pnrep[m] = 1/(1+exp(-etarep[m]));
    // pnrep[m] = beta_rng(an, bn);
    Nrep[m] = binomial_rng(Nt, pnrep[m]);
    // log_lik[m] = binomial_lpmf(rep_vector(N[m],1) | Nt, pnrep[m]);
  }
  for (m in 1:M) {
    log_lik[m] = binomial_lpmf(NP3[m] | Nt, pn3[m]);
  }
}





