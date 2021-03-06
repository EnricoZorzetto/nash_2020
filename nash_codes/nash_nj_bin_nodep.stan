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
 real<lower=0, upper=1> pn;
}
model {
  // Priors
  pn ~ beta(4, 4);
  // ma2 ~ beta(2, 2);
  // bn ~ gamma(10,1); // I want something fairly uninformative
  // sx ~ normal(0, 2);
  // sy ~ normal(0, 1);
  for (m in 1:M) {
    // pn[m] ~ beta(an[m], bn);  
  target += binomial_lpmf(N[m] | Nt, pn);
  }
}
generated quantities {
  vector[M] log_lik;
  vector[Mgen] Nrep;
  // vector[M] pnrep;
  // vector[Nt] predpmf;
  for (m in 1:Mgen) {
    Nrep[m] = binomial_rng(Nt, pn);
    // log_lik[m] = binomial_lpmf(rep_vector(N[m],1) | Nt, pnrep[m]);
  }  
  for (m in 1:M) {
    log_lik[m] = binomial_lpmf(N[m] | Nt, pn);
  }
}





