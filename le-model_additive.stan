parameters {

  // matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age;
  vector[nk_incrate_time] coef_incrate_time;
  vector[nk_incrate_age-1] param_incrate_age;
  vector[nk_natmx_time] coef_natmx_time;
  vector[nk_natmx_age-1] param_natmx_age;
  vector<upper=0>[STEPS_time-artstart_tIDX] dt_log_artrr;
  // real<lower=0> sigma_incrate_time_age;
  real<lower=0> sigma_incrate_time;
  real<lower=0> sigma_incrate_age;
  real<lower=0> sigma_natmx_time;
  real<lower=0> sigma_natmx_age;
  real<lower=0> sigma_art;
}
transformed parameters{

  vector[nk_incrate_age] coef_incrate_age;
  vector[nk_natmx_age] coef_natmx_age;
  matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age;

  for(i in 1:nk_natmx_age)
    if (i < fixcoef_natmx_age){
      coef_incrate_age[i] <- param_incrate_age[i];
      coef_natmx_age[i] <- param_natmx_age[i];
    } else if (i == fixcoef_natmx_age) {
      coef_incrate_age[i] <- -sum(param_incrate_age);
      coef_natmx_age[i] <- -sum(param_natmx_age);
    } else {
      coef_incrate_age[i] <- param_incrate_age[i-1];
      coef_natmx_age[i] <- param_natmx_age[i-1];
    }

  for(i in 1:nk_incrate_time)
    for(j in 1:nk_incrate_age)
      coef_incrate_time_age[i, j] <- coef_incrate_time[i] + coef_incrate_age[j];
}
model {

  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  // sigma_incrate_time_age ~ cauchy(0, 2.5);
  sigma_incrate_time ~ cauchy(0, 2.5);
  sigma_incrate_age ~ cauchy(0, 2.5);
  sigma_natmx_time ~ cauchy(0, 2.5);
  sigma_natmx_age ~ cauchy(0, 2.5);
  sigma_art ~ cauchy(0, 2.5);

  //////////////////////
  //  Spline penalty  //
  //////////////////////

  {

    D_incrate_time * coef_incrate_time ~ normal(0, sigma_incrate_time);
    D_incrate_age * coef_incrate_age ~ normal(0, sigma_incrate_age);
    
    D_natmx_time * coef_natmx_time ~ normal(0, sigma_natmx_time);
    D_natmx_age * coef_natmx_age ~ normal(0, sigma_natmx_age);
    D_art * dt_log_artrr ~ normal(0, sigma_art);
  }
