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

  // incidence model
  int<lower=5> nk_incrate_time;
  int<lower=5> nk_incrate_age;

  matrix[STEPS_time, nk_incrate_time] X_incrate_time;
  matrix[STEPS_time-1, nk_incrate_time] Xmid_incrate_time;
  matrix[STEPS_age, nk_incrate_age] X_incrate_age;
  matrix[STEPS_age-1, nk_incrate_age] Xmid_incrate_age;

  int<lower=0> pen_ord_incrate;

  matrix[nk_incrate_time-pen_ord_incrate, nk_incrate_time] D_incrate_time;
  matrix[nk_incrate_age-pen_ord_incrate, nk_incrate_age] D_incrate_age;
  matrix[nk_incrate_time*nk_incrate_age, nk_incrate_time*nk_incrate_age] Pcar_prec_incrate;


  // non-HIV mortality model
  int<lower=5> nk_natmx_time;
  int<lower=5> nk_natmx_age;

  matrix[STEPS_time, nk_natmx_time] X_natmx_time;
  matrix[STEPS_time-1, nk_natmx_time] Xmid_natmx_time;
  matrix[STEPS_age, nk_natmx_age] X_natmx_age;
  matrix[STEPS_age-1, nk_natmx_age] Xmid_natmx_age;

  int<lower=0> pen_ord_natmx_time;
  int<lower=0> pen_ord_natmx_age;

  matrix[nk_natmx_time-pen_ord_natmx_time, nk_natmx_time] D_natmx_time;
  matrix[nk_natmx_age-pen_ord_natmx_age, nk_natmx_age] D_natmx_age;

  int<lower=1, upper=nk_natmx_time> fixcoef_natmx_time;
  int<lower=1, upper=nk_natmx_age> fixcoef_natmx_age;


  // ART model
  
  // matrix[STEPS_time, nk_art] X_art;
  // matrix[STEPS_time-1, nk_art] Xmid_art;

  int<lower=0> pen_ord_art;
  matrix[STEPS_time-artstart_tIDX-pen_ord_art, STEPS_time-artstart_tIDX] D_art;
  
  
  // HIV survival model
  matrix[STEPS_time-1, STEPS_age-1] hivmx_dur_a0;     // sequenced [1:DUR, 1:STEPS_age]
  matrix[STEPS_time-1, STEPS_age-1] hivsurv_dur_a0;   // sequenced [1:DUR, 1:STEPS_age]
  matrix[STEPS_time-1, STEPS_age-1] hivmxMID_dur_a0;   // sequenced [1:(STEPS_time-1), 1:STEPS_age]

}
transformed data {
}
