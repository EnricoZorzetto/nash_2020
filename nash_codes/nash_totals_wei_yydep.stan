
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  real<lower=0> totals[M]; // seasonal totals
  real XN[M];  // latitude of the pressure ridge
  real YN[M];  // longitude of the pressure ridge
  real XNgen[Mgen];  // latitude of the pressure ridge
  real YNgen[Mgen];  // longitude of the pressure ridge
}
transformed data{
  int<lower=0, upper=1> are_totals_nonzero[M]; 
  int<lower=0> Mwets;
  Mwets = 0;
  for (m in 1:M){
    if (totals[m] >1e-6){
      are_totals_nonzero[m] = 1;
      Mwets = Mwets + 1;
    } else {
      are_totals_nonzero[m] = 0;
    }
  }
}
parameters {
  real<lower=0, upper=1> p0;
  real<lower=0> w;
  real c0;
  real cy;
}
transformed parameters{
  real<lower=0> C[M];
  for(m in 1:M){
    C[m] = exp(c0 + cy*YN[m]);
  }
}
model {
  // Priors
  p0 ~ beta(5,100);
  c0 ~ normal(3.5, 1);
  cy ~ normal(0, 0.1);
  for (m in 1:M) {
    if (totals[m] > 1e-6){
      totals[m] ~ weibull( w, C[m]);
      are_totals_nonzero[m]~bernoulli(p0);
    } else {
      are_totals_nonzero[m]~bernoulli(p0);
    }
  }
}
generated quantities {
  vector[Mwets] log_lik;
  // vector[Mgen] totalsrep;
  int<lower=1> count;
  vector[Mgen] Crep;
  for (m in 1:Mgen) {
    Crep[m] = exp(c0 + cy*YNgen[m]);
    // totalsrep[m] = weibull_rng(w, Crep[m]);
    // log_lik[m] = normal_lpdf(totals[m] | mu[m], sigma);
  }
  count = 1;
  for (m in 1:M){
    if (totals[m]>1e-6){
      log_lik[count] = weibull_lpdf(totals[m] | w, C[m]);
      count = count + 1;
    }
  }
}




