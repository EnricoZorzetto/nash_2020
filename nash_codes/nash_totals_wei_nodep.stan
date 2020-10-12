
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
  real<lower=0> C;
}
model {
  // Priors
  p0 ~ beta(5,100);
  for (m in 1:M) {
    if (totals[m] > 1e-6){
      totals[m] ~ weibull( w, C);
      are_totals_nonzero[m]~bernoulli(p0);
    } else {
      are_totals_nonzero[m]~bernoulli(p0);
    }
  }
}
generated quantities {
  vector[Mwets] log_lik;
  // vector[Mgen] totalsrep;
  // vector[Mgen] Crep;
  int<lower=1> count;
  count = 1;
  for (m in 1:M) {
    // Crep[m] = exp(mc0 + mcy*YN[m]);
    // totalsrep[m] = weibull_rng(w, Crep[m]);
    // log_lik[m] = normal_lpdf(totals[m] | mu[m], sigma);
    if (totals[m]>1e-6){
    log_lik[count] = weibull_lpdf(totals[m] | w, C);
    count = count + 1;
    }
  }
}




