  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years
  int<lower = 1> Nt; // number of years
  int<lower=0> N[M];  // number of events - sample
  real XN[M];  // number of events - sample
  real YN[M];  // number of events - sample
  real XNgen[Mgen];  // number of events - sample
  real YNgen[Mgen];  // number of events - sample
}
parameters {
  // real sx;
  real n0;
  real ny;
}
transformed parameters{
 // real<lower=0> an[M];
 // real<lower=0> bn[M];
 real<lower=0, upper=1> pn[M];
 // real ma;
 real eta[M];
 // real mnv[M];
 // ma = 0.5*ma2;
 for (m in 1:M){
   eta[m] = n0 + ny*YN[m];
   pn[m] = 1/(1+exp(-eta[m]));
 }
}
model {
  // Priors
  // mn0 ~ beta(2, 2);
  // ma2 ~ beta(2, 2);
  // bn ~ gamma(10,1); // I want something fairly uninformative
  // sx ~ normal(0, 2);
  // sy ~ normal(0, 1);
  n0 ~ normal(0, 1);
  ny ~ normal(0, 0.5);
  for (m in 1:M) {
    // pn[m] ~ beta(an[m], bn);  
  target += binomial_lpmf(N[m] | Nt, pn[m]);
  }
}
generated quantities {
  vector[M] log_lik;
  vector[Mgen] Nrep;
  vector[Mgen] pnrep;
  vector[Mgen] etarep;
  // vector[Nt] predpmf;
  for (m in 1:Mgen) {
    etarep[m] = n0 + ny*YNgen[m];
    pnrep[m] = 1/(1+exp(-etarep[m]));
    // pnrep[m] = beta_rng(an, bn);
    Nrep[m] = binomial_rng(Nt, pnrep[m]);
    // log_lik[m] = binomial_lpmf(rep_vector(N[m],1) | Nt, pnrep[m]);
  }
  for (m in 1:M) {
    log_lik[m] = binomial_lpmf(N[m] | Nt, pn[m]);
  }
}




