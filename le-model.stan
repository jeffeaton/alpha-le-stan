parameters {

  matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age;
  vector[nk_natmx_time] coef_natmx_time;
  vector[nk_natmx_age-1] param_natmx_age;
  vector<upper=0>[STEPS_time-artstart_tIDX] dt_log_artrr;
  real<lower=0> sigma_incrate_time_age;
  real<lower=0> sigma_incrate_time;
  real<lower=0> sigma_incrate_age;
  real<lower=0> sigma_natmx_time;
  real<lower=0> sigma_natmx_age;
  real<lower=0> sigma_art;
}
transformed parameters{

  vector[nk_incrate_time] coef_incrate_time;
  row_vector[nk_incrate_age] coef_incrate_age;
  vector[nk_natmx_age] coef_natmx_age;

  coef_incrate_time <- row_means(coef_incrate_time_age);
  coef_incrate_age <- col_means(coef_incrate_time_age);
    
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

  matrix[STEPS_time-1, STEPS_age-1] incrateMID_time_age;
  matrix[STEPS_time, STEPS_age] cumavoid_time_age;
  matrix[STEPS_time-1, STEPS_age-1] cumavoidMID_time_age;
  matrix[STEPS_time, STEPS_age] natmx_time_age;
  matrix[STEPS_time, STEPS_age] natsurv_time_age;
  vector[STEPS_time] artrr;
  vector[STEPS_time-1] artrr_MID;

  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  sigma_incrate_time_age ~ cauchy(0, 2.5);
  sigma_incrate_time ~ cauchy(0, 2.5);
  sigma_incrate_age ~ cauchy(0, 2.5);
  sigma_natmx_time ~ cauchy(0, 2.5);
  sigma_natmx_age ~ cauchy(0, 2.5);
  sigma_art ~ cauchy(0, 2.5);

  //////////////////////
  //  Spline penalty  //
  //////////////////////

  {
    matrix[nk_incrate_time, nk_incrate_age] resid_incrate_time_age;
    vector[nk_incrate_time*nk_incrate_age] vec_resid_incrate_time_age;

    for(j in 1:nk_incrate_age)
      for(i in 1:nk_incrate_time)
	resid_incrate_time_age[i,j] <- coef_incrate_time_age[i,j] - coef_incrate_time[i] - coef_incrate_age[j];

    vec_resid_incrate_time_age <- to_vector(resid_incrate_time_age);
    increment_log_prob(-nk_incrate_time*nk_incrate_age*log(sigma_incrate_time_age) -
		       1/(2*sigma_incrate_time_age*sigma_incrate_time_age) * (vec_resid_incrate_time_age' * Pcar_prec_incrate * vec_resid_incrate_time_age));
    
    D_natmx_time * coef_natmx_time ~ normal(0, sigma_natmx_time);
    D_natmx_age * coef_natmx_age ~ normal(0, sigma_natmx_age);
    D_art * dt_log_artrr ~ normal(0, sigma_art);
  }
  
  ///////////////////////////////////////////////////////
  //  Construct incidence rate and mortality matrices  //
  ///////////////////////////////////////////////////////

  incrateMID_time_age <- exp(Xmid_incrate_time * coef_incrate_time_age * Xmid_incrate_age');
  cumavoid_time_age <- exp(-dt*diagCumSum(incrateMID_time_age));
  cumavoidMID_time_age <- block(cumavoid_time_age, 1, 1, STEPS_time-1, STEPS_age-1) .* exp(-dt/2*incrateMID_time_age);

  natmx_time_age <- exp(X_natmx_time * coef_natmx_time) * exp(X_natmx_age * coef_natmx_age)';
  natsurv_time_age <- exp(-dt*diagCumSum(exp(Xmid_natmx_time * coef_natmx_time) * exp(Xmid_natmx_age * coef_natmx_age)'));

  {
    vector[STEPS_time-artstart_tIDX] log_artrr;
    log_artrr <- dt*cumulative_sum(dt_log_artrr);

    for(i in 1:(STEPS_time-artstart_tIDX)){
      artrr[artstart_tIDX+i] <-  exp(log_artrr[i]);
      artrr_MID[artstart_tIDX+i-1] <- exp(log_artrr[i] - dt/2*dt_log_artrr[i]);
    }
  }

  ///////////////////////////////////////
  //  Calculate individual likelihood  //
  ///////////////////////////////////////

  increment_log_prob(calc_ll_cohexit(coh_cIDX, coh_minexpose_tIDX, coh_maxexpose_tIDX, coh_nexit,
  				     exdat_tIDX, exdat_minexpose_tIDX, exdat_maxexpose_tIDX, exdat_ndat,
  				     aggr_exposestart_tIDX, aggr_exposeend_tIDX, aggr_death, aggr_hivpos, aggr_nrepl,
  				     cumavoid_time_age, cumavoidMID_time_age, incrateMID_time_age,
  				     hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
  				     artrr, artrr_MID, artstart_tIDX,
  				     natsurv_time_age, natmx_time_age, dt));
}
