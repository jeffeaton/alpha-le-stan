data {

  // STATE SPACE PARAMETERS
  real<lower=0> dt;
  int<lower=1> STEPS_time;
  int<lower=1> STEPS_age;
  int<lower=1, upper=STEPS_time> artstart_tIDX;


  // COHORT DATA
  int<lower=1> NCOH;

  int coh_cIDX[NCOH]; 
  int<lower=1, upper=STEPS_time> coh_minexpose_tIDX[NCOH];
  int<lower=1, upper=STEPS_time-1> coh_maxexpose_tIDX[NCOH];
  int<lower=1> coh_nexit[NCOH];

    // EXIT DATA
  int<lower=NCOH> NEXIT;

  int<lower=1, upper=STEPS_time> exdat_tIDX[NEXIT];
  int<lower=1, upper=STEPS_time> exdat_minexpose_tIDX[NEXIT];
  int<lower=1, upper=STEPS_time> exdat_maxexpose_tIDX[NEXIT];
  int<lower=1> exdat_ndat[NEXIT];
  
  // AGGREGATE INDIVIDUAL DATA
  int<lower=NEXIT> NAGGR;

  int<lower=1, upper=STEPS_time> aggr_exposestart_tIDX[NAGGR];
  int<lower=1, upper=STEPS_time> aggr_exposeend_tIDX[NAGGR];
  int<lower=0, upper=1> aggr_hivpos[NAGGR];
  int<lower=0, upper=1> aggr_death[NAGGR];
  vector[NAGGR] aggr_nrepl;


  // MODEL PARAMETERS
  int<lower=5> nk_time;
  int<lower=5> nk_age;
  int<lower=5> nk_art;
  int<lower=5> nk_natmx;

  int<lower=1, upper=nk_time> fixcoef_time_idx;
  int<lower=1, upper=nk_age> fixcoef_age_idx;

  matrix[STEPS_time, nk_time] X_time;
  matrix[STEPS_age, nk_age] X_age;
  matrix[STEPS_time, nk_art] X_art;
  matrix[STEPS_time, nk_natmx] X_natmx;
  matrix[STEPS_time-1, nk_time] Xmid_time;
  matrix[STEPS_age-1, nk_age] Xmid_age;
  matrix[STEPS_time-1, nk_art] Xmid_art;
  matrix[STEPS_time-1, nk_natmx] Xmid_natmx;

  // matrix[nk_time-1, nk_time] P_time;
  matrix[nk_age-1, nk_age] P_age;
  // matrix[nk_art-1, nk_art] P_art;
  matrix[STEPS_time-artstart_tIDX-1, STEPS_time-artstart_tIDX] P_art;
  matrix[nk_natmx-1, nk_natmx] P_natmx;

  matrix[nk_time*nk_age, nk_time*nk_age] Pcar_prec_incrate;

  // vector[nk_time] coef_incrate_time;
  // vector[nk_age] coef_incrate_age;
  // vector[nk_time] coef_natmx_time;
  // vector[nk_age] coef_natmx_age;
  // vector[nk_art] coef_art;

  // real<lower=0> sigma2_incrate_time;
  // real<lower=0> sigma2_incrate_age;
  // real<lower=0> sigma2_incrate_time_age;
  // real<lower=0> sigma2_natmx_time;
  // real<lower=0> sigma2_natmx_age;
  // real<lower=0> sigma2_art;

  matrix[STEPS_time-1, STEPS_age-1] hivmx_dur_a0;     // sequenced [1:DUR, 1:STEPS_age]
  matrix[STEPS_time-1, STEPS_age-1] hivsurv_dur_a0;   // sequenced [1:DUR, 1:STEPS_age]
  matrix[STEPS_time-1, STEPS_age-1] hivmxMID_dur_a0;   // sequenced [1:(STEPS_time-1), 1:STEPS_age]

}
transformed data {
}
