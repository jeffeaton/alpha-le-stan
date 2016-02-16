parameters {

  matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age;
  vector[nk_natmx_time] coef_natmx_time;
  vector[nk_natmx_age-1] param_natmx_age;
  vector<upper=0>[STEPS_time-artstart_tIDX] dt_log_artrr;
  real<lower=0> sigma_incrate_time_age;
  real<lower=0> sigma_natmx_time;
  real<lower=0> sigma_natmx_age;
  real<lower=0> sigma_art;
}
transformed parameters{

  vector[nk_natmx_age] coef_natmx_age;

  for(i in 1:nk_natmx_age)
    if (i < fixcoef_natmx_age){
      coef_natmx_age[i] <- param_natmx_age[i];
    } else if (i == fixcoef_natmx_age) {
      coef_natmx_age[i] <- -sum(param_natmx_age);
    } else {
      coef_natmx_age[i] <- param_natmx_age[i-1];
    }
}
model {


  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  sigma_incrate_time_age ~ cauchy(0, 2.5);
  sigma_natmx_time ~ cauchy(0, 2.5);
  sigma_natmx_age ~ cauchy(0, 2.5);
  sigma_art ~ cauchy(0, 2.5);

  //////////////////////
  //  Spline penalty  //
  //////////////////////

  {
    vector[nk_incrate_time*nk_incrate_age] vec_coef_incrate_time_age;

    vec_coef_incrate_time_age <- to_vector(coef_incrate_time_age);
    increment_log_prob(-nk_incrate_time*nk_incrate_age*log(sigma_incrate_time_age) -
		       1/(2*sigma_incrate_time_age*sigma_incrate_time_age) * (vec_coef_incrate_time_age' * Pcar_prec_incrate * vec_coef_incrate_time_age));
    
    D_natmx_time * coef_natmx_time ~ normal(0, sigma_natmx_time);
    D_natmx_age * coef_natmx_age ~ normal(0, sigma_natmx_age);
    D_art * dt_log_artrr ~ normal(0, sigma_art);
  }
  
