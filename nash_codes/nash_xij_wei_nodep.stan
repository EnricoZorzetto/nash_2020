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
}
transformed data{
  int<lower=0> nrepobs;
  nrepobs = Nt*Mgen;
}
parameters {
  real<lower=0> w;
  real<lower=0> C;
  real<lower=0, upper=1> pn;
}
model {
  pn ~ beta(4, 4);
  w ~ gamma(0.7*10, 10);
  C ~ gamma(7*2, 2);
    for (m in 1:M) {
    target += binomial_lpmf(N[m] | Nt, pn);
      for (j in 1:Nt) {
        if (y[m,j] >1e-6) {
           y[m,j] ~ weibull( w, C);
        }
      }
    }
}
generated quantities {
  vector[M] log_lik_nj;
  int<lower=0> Nrep[Mgen];
  vector[totnwets] log_lik; // log like for the N model only
  vector[Mgen] Crep;
  vector[Mgen] wrep;
  vector[nrepobs] xijrep;
  vector[Mgen] maxrep;
  int<lower=0> count;
  int<lower=0> startm;
  int<lower=0> endm;
  
  for (m in 1:Mgen) {
    Crep[m] = C;
    wrep[m] = w;
    Nrep[m] = binomial_rng(Nt, pn);
    for (j in 1:Nt) {
      if ( j <= Nrep[m]){
        xijrep[ (m-1)*Nt + j ] = weibull_rng(wrep[m], Crep[m]);
      } else {
        xijrep[ (m-1)*Nt + j ] = 0;
      }
    }
    startm = (m-1)*Nt+1;
    endm = (m-1)*Nt+Nrep[m];
    maxrep[m] = max(xijrep[startm:endm]);
  }  
  count = 0;
  for (m in 1:M) {
    log_lik_nj[m] = binomial_lpmf(N[m] | Nt, pn);
    for (j in 1:Nt) {
      if (y[m,j] > 1e-6){
        count = count + 1;
        log_lik[count] = weibull_lpdf(y[m,j] | w, C);
        }
  }
}
}





