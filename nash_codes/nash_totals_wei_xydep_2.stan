  // Hierarchical model for extreme events
  functions {
    // compound Weibull with atom of probability in zero
  // real p0weibull_lpdf(vector y, real C, real w, real p0) {
  //   // generalised Pareto log pdf 
  //   real mylpdf;
  //   int N = rows(y);
  //   for (i in 1:N){
  //   if (fabs(y[i]) > 1e-15){
  //     mylpdf = mylpdf + log(1-p0) + weibull_lpdf(y | w, C);
  //   } else {
  //     mylpdf = mylpdf + log(p0);
  //   }
  //   }
  //   return mylpdf;
  // }
    real p0weibull_lpdf(real y, real w, real C, real p0) {
    // generalised Pareto log pdf 
    real mylpdf;
    // int N = rows(y);
    // for (i in 1:N){
    if (fabs(y) > 1e-15){
      mylpdf = log(1-p0) + weibull_lpdf(y | w, C);
    } else {
      mylpdf = log(p0);
    }
    // }
    return mylpdf;
  }
  }
//   real p0weibull_rng(real w, real C, real p0) {
//     // generalised Pareto rng
//     real Fi;
//     Fi = uniform_rng(0,1);
//     if (Fi < p0){
//       return 0
//     } else {
//       return p0 + (1-p0)*weibull()
//     }
//       reject("sigma<=0; found sigma =", sigma)
//     if (fabs(k) > 1e-15)
//       return ymin + (uniform_rng(0,1)^-k -1) * sigma / k;
//     else
//       return ymin - sigma*log(uniform_rng(0,1)); // limit k->0
//   } 
// }

data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  real<lower=0> totals[M]; // seasonal totals
  real XN[M];  // latitude of the pressure ridge
  real YN[M];  // longitude of the pressure ridge
  real XNgen[Mgen];  // latitude of the pressure ridge
  real YNgen[Mgen];  // longitude of the pressure ridge
}
parameters {
  real<lower=0, upper=1> p0;
  real<lower=0> w;
  real c0;
  real cy;
  real cx;
  // real muy;
  // real sx;
}
// transformed parameters{
//  real<lower=0> mu[M];
//  for(m in 1:M){
//  mu[m] = mu0 + muy*YN[m]*1/(1+exp(-sx*XN[m]));
//  }
// }

transformed parameters{
 // real<lower=0> mc[M];
 real<lower=0> C[M];
 // real eta[M];

 for(m in 1:M){
 // mc[m] = mc0 + mcy*YN[m]*1/(1+exp(-sx*XN[m]));
 C[m] = exp(c0 + cy*YN[m]*1/(1+exp(-cx*XN[m])));
 }
}
model {
  // Priorss
  p0 ~ beta(5,100);
  // c0 ~ normal(0, 1);
  c0 ~ normal(3.5, 1);
  cy ~ normal(0, 0.5);
  cx ~ normal(0, 1);
    for (m in 1:M) {
    // totals[m] ~ normal( mu[m], sigma);
   totals[m] ~ p0weibull( w, C[m], p0);
   // totals[m] ~ p0weibull( w, C[m], p0);
    }
}
generated quantities {
  vector[M] log_lik;
  // vector[Mgen] totalsrep;
  vector[Mgen] Crep;
  for (m in 1:Mgen){
    Crep[m] = exp(c0 + cy*YNgen[m]*1/(1+exp(-cx*XNgen[m])));
    // totalsrep[m] = weibull_rng(w, Crep[m]);
    // log_lik[m] = normal_lpdf(totals[m] | mu[m], sigma);
  }
  for (m in 1:M){
    log_lik[m] = p0weibull_lpdf(totals[m] | w, C[m], p0);
  }
}




