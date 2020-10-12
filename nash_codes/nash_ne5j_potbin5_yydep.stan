  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years
  int<lower = 1> Nt; // number of years
  int<lower=0> NP5[M];  // number of events - sample
  real XN[M];  // number of events - sample
  real YN[M];  // number of events - sample
  real XNgen[Mgen];  // number of events - sample
  real YNgen[Mgen];  // number of events - sample
}
parameters {
  // real sx;
  real n50;
  real n5y;
}
transformed parameters{
 // real<lower=0> an[M];
 // real<lower=0> bn[M];
 real<lower=0, upper=1> pn5[M];
 // real ma;
 real eta5[M];
 // real mnv[M];
 // ma = 0.5*ma2;
 for (m in 1:M){
   eta5[m] = n50 + n5y*YN[m];
   pn5[m] = 1/(1+exp(-eta5[m]));
 }
}
model {
  // Priors
  // mn0 ~ beta(2, 2);
  // ma2 ~ beta(2, 2);
  // bn ~ gamma(10,1); // I want something fairly uninformative
  // sx ~ normal(0, 2);
  // sy ~ normal(0, 1);
  for (m in 1:M) {
    // pn[m] ~ beta(an[m], bn);  
  target += binomial_lpmf(NP5[m] | Nt, pn5[m]);
  }
}
generated quantities {
  vector[M] log_lik;
  vector[Mgen] Nrep;
  vector[Mgen] pnrep;
  vector[Mgen] etarep;
  // vector[Nt] predpmf;
  for (m in 1:Mgen) {
    etarep[m] = n50 + n5y*YNgen[m];
    pnrep[m] = 1/(1+exp(-etarep[m]));
    // pnrep[m] = beta_rng(an, bn);
    Nrep[m] = binomial_rng(Nt, pnrep[m]);
    // log_lik[m] = binomial_lpmf(rep_vector(N[m],1) | Nt, pnrep[m]);
  }
  for (m in 1:M) {
    log_lik[m] = binomial_lpmf(NP5[m] | Nt, pn5[m]);
  }
}




