functions {

  matrix diagCumSum(matrix x){

    matrix[rows(x)+1, cols(x)+1] val;

    for(i in 1:rows(val))
      val[i,1] <- 0;
    for(j in 2:cols(val)){
      val[1,j] <- 0;
      for(i in 2:rows(val))
        val[i,j] <- val[i-1,j-1] + x[i-1,j-1];
    }

    return val;
  }

  // real hivsurv(real tao, real a0){
  //   real k;
  //   real Lambda;

  //   k <- 2.3;
  //   Lambda <- 166.0/(k*a0^0.53);
  //   return exp(-(tao/Lambda)^k);
  // }

  // vector hivsurv(vector tao, vector a0){
  //   real k;
  //   vector[rows(tao)] Lambda;

  //   k <- 2.3;
  //   Lambda <- 166.0 ./(k*exp(0.53*log(a0)));
  //   return exp(-exp(k*log(tao ./Lambda)));
  // }

  // real hivsurvhaz(real tao, real a0){
  //   real k;
  //   real Lambda;

  //   k <- 2.3;
  //   Lambda <- 166.0/(k*exp(0.53*log(a0)));

  //   return k/Lambda*(tao/Lambda)^(k-1);
  // }

  // vector hivsurvhaz(vector tao, vector a0){
  //   real k;
  //   vector[rows(tao)] Lambda;

  //   k <- 2.3;
  //   Lambda <- 166.0 ./(k*exp(0.53*log(a0)));
  //   return k ./ Lambda .* exp((k-1)*log(tao ./ Lambda));
  // }

  vector calc_log_phivsurv(int tIDX, int aIDX, int taoIDX, int exposeDUR, int exit_tIDX,
                           matrix log_hivsurv_REVdur_a0, matrix log_hivmxMID_dur_a0, vector log_artrr_MID, int artstart_tIDX, real logdt){
    vector[exposeDUR] log_phivsurv;

    if(exit_tIDX > artstart_tIDX){

      int i_tIDX;
      int i_aIDX;

      log_phivsurv <- rep_vector(0, exposeDUR);

      for(ii in 1:exposeDUR){
        i_tIDX <- tIDX + ii - 1;
        i_aIDX <- aIDX + ii - 1;

        log_phivsurv[ii] <- -exp(log_sum_exp(sub_col(log_hivmxMID_dur_a0, 1, i_aIDX, exit_tIDX - i_tIDX) + segment(log_artrr_MID, i_tIDX, exit_tIDX - i_tIDX)) + logdt);
      }
    } else
      log_phivsurv <- diagonal(block(log_hivsurv_REVdur_a0, taoIDX, aIDX, exposeDUR, exposeDUR));

    return(log_phivsurv);
  }

  vector calc_log_hivmx(int tIDX, int aIDX, int taoIDX, int exposeDUR, int exit_tIDX,
                        matrix log_hivmx_REVdur_a0, vector log_artrr, int artstart_tIDX){
    vector[exposeDUR] log_hivmx;
    log_hivmx <- diagonal(block(log_hivmx_REVdur_a0, taoIDX, aIDX, exposeDUR, exposeDUR));
    if(tIDX + exposeDUR > artstart_tIDX)
      log_hivmx <- log_hivmx + log_artrr[exit_tIDX];
    return(log_hivmx);
  }
}
data {

  // STATE SPACE PARAMETERS
  real<lower=0> dt;
  int<lower=1> STEPS_time;
  int<lower=1> STEPS_age;
  int<lower=1, upper=STEPS_time> artstart_tIDX;


  // INDIVIDUAL DATA
  int<lower=1> N;

  int<lower=1, upper=STEPS_time> entry_tIDX[N];
  int<lower=1, upper=STEPS_time> exit_tIDX[N];
  int<lower=1, upper=STEPS_time> exposestart_tIDX[N];

  int<lower=1, upper=STEPS_age> entry_aIDX[N];
  int<lower=1, upper=STEPS_age> exit_aIDX[N];
  int<lower=1, upper=STEPS_age> exposestart_aIDX[N];

  int<lower=0> expose_DUR[N];


  int<lower=0, upper=1> hivpos[N];
  int<lower=0, upper=1> death[N];


  // MODEL PARAMETERS
  int<lower=5> nk_time;
  int<lower=5> nk_age;
  int<lower=5> nk_art;

  int<lower=1, upper=nk_age> fixcoef_age_idx;

  matrix[STEPS_time, nk_time] X_time;
  matrix[STEPS_age, nk_age] X_age;
  matrix[STEPS_time, nk_art] X_art;
  matrix[STEPS_time-1, nk_time] Xmid_time;
  matrix[STEPS_age-1, nk_age] Xmid_age;
  matrix[STEPS_time-1, nk_art] Xmid_art;

  matrix[nk_time-1, nk_time] P_time;
  matrix[nk_age-1, nk_age] P_age;
  matrix[nk_art-1, nk_art] P_art;

  // vector[nk_time] coef_incrate_time;
  // vector[nk_age] coef_incrate_age;
  // vector[nk_time] coef_natmx_time;
  // vector[nk_age] coef_natmx_age;
  // vector[nk_art] coef_art;

  // real<lower=0> sigma2_art;

  matrix[STEPS_time, STEPS_age] log_hivmx_REVdur_a0;     // sequenced [STEPS_time:1, 1:STEPS_age]
  matrix[STEPS_time, STEPS_age] log_hivsurv_REVdur_a0;   // sequenced [STEPS_time:1, 1:STEPS_age]
  matrix[STEPS_time-1, STEPS_age] log_hivmxMID_dur_a0;   // sequenced [1:(STEPS_time-1), 1:STEPS_age]  

}
transformed data{

  real logdt;

  int<lower=1, upper=STEPS_time+1> exposestart_taoIDX[N];  // STEPS_time+1 if expose_DUR=0
  
  int<lower=0, upper=STEPS_time> lefttrunc_expose_tIDX[N];
  int<lower=0, upper=STEPS_age> lefttrunc_expose_aIDX[N];
  int<lower=0, upper=STEPS_age> lefttrunc_expose_taoIDX[N];
  int<lower=0, upper=min(STEPS_time, STEPS_age)> lefttrunc_expose_DUR[N];

  for(i in 1:N){
    exposestart_taoIDX[i] <- exposestart_tIDX[i] + STEPS_time - exit_tIDX[i] + 1;

    lefttrunc_expose_DUR[i] <- min(entry_tIDX[i], entry_aIDX[i]) - 1;
    lefttrunc_expose_tIDX[i] <- entry_tIDX[i] - lefttrunc_expose_DUR[i];
    lefttrunc_expose_aIDX[i] <- entry_aIDX[i] - lefttrunc_expose_DUR[i];
    lefttrunc_expose_taoIDX[i] <- lefttrunc_expose_tIDX[i] + STEPS_time - entry_tIDX[i] + 1;
  }

  logdt <- log(dt);
}
parameters {
  real<lower=0> sigma2_incrate_time;
  real<lower=0> sigma2_incrate_age;
  real<lower=0> sigma2_natmx_time;
  real<lower=0> sigma2_natmx_age;
  real<lower=0> sigma2_art;

  vector[nk_time] coef_incrate_time;
  vector[nk_age-1] param_incrate_age;
  vector[nk_time] coef_natmx_time;
  vector[nk_age-1] param_natmx_age;
  vector[nk_art] coef_art;
  // real<lower=-25, upper=25> param_art;

}
transformed parameters{
  vector[nk_age] coef_incrate_age;
  vector[nk_age] coef_natmx_age;

//  vector[nk_art] coef_art;

// coef_art <- rep_vector(param_art, nk_art);

  for(i in 1:nk_age)
    if (i < fixcoef_age_idx){
      coef_incrate_age[i] <- param_incrate_age[i];
      coef_natmx_age[i] <- param_natmx_age[i];
    } else if (i == fixcoef_age_idx) {
      coef_incrate_age[i] <- 0.0;
      coef_natmx_age[i] <- 0.0;      
    } else {
      coef_incrate_age[i] <- param_incrate_age[i-1];
      coef_natmx_age[i] <- param_natmx_age[i-1];
    }
}
model {

  matrix[STEPS_time, STEPS_age] log_incrate_time_age;
  matrix[STEPS_time, STEPS_age] log_cumavoid_time_age;
  matrix[STEPS_time, STEPS_age] log_natmx_time_age;
  matrix[STEPS_time, STEPS_age] log_natsurv_time_age;
  vector[STEPS_time] log_artrr;
  vector[STEPS_time-1] log_artrr_MID;

  real log_psurventry;
  real log_phivn;
  real log_phivp;

  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  1/sigma2_incrate_time ~ gamma(1.0, 0.0005);
  1/sigma2_incrate_age ~ gamma(1.0, 0.0005);
  1/sigma2_natmx_time ~ gamma(1.0, 0.0005);
  1/sigma2_natmx_age ~ gamma(1.0, 0.0005);
  1/sigma2_art ~ gamma(1.0, 0.0005);

  //////////////////////
  //  Spline penalty  //
  //////////////////////
  
  P_time * coef_incrate_time ~ normal(0, sqrt(sigma2_incrate_time));
  P_age * coef_incrate_age ~ normal(0, sqrt(sigma2_incrate_age));
  P_time * coef_natmx_time ~ normal(0, sqrt(sigma2_natmx_time));
  P_age * coef_natmx_age ~ normal(0, sqrt(sigma2_natmx_age));
  P_art * coef_art ~ normal(0, sqrt(sigma2_art));

  ///////////////////////////////////////////////////////
  //  Construct incidence rate and mortality matrices  //
  ///////////////////////////////////////////////////////

  log_incrate_time_age <- log(exp(X_time * coef_incrate_time) * exp(X_age * coef_incrate_age)');
  log_cumavoid_time_age <- -dt*diagCumSum(exp(Xmid_time * coef_incrate_time) * exp(Xmid_age * coef_incrate_age)');

  log_natmx_time_age <- log(exp(X_time * coef_natmx_time) * exp(X_age * coef_natmx_age)');
  log_natsurv_time_age <- -dt*diagCumSum(exp(Xmid_time * coef_natmx_time) * exp(Xmid_age * coef_natmx_age)');

  log_artrr <- X_art * coef_art;
  log_artrr_MID <- Xmid_art * coef_art;

// print("coef_art ", coef_art);
// print("log_artrr_REV ", log_artrr_REV);

// print(block(log_incrate_time_age, 1, 1, 5, 5));
// print("coef_incrate_time ", coef_incrate_time);
// print("coef_incrate_age ", coef_incrate_age);
// print("coef_natmx_time ", coef_natmx_time);
// print("coef_natmx_age ", coef_natmx_age);
// print("coef_art ", coef_art, ", sigma2_art ", sigma2_art);


  ///////////////////////////////////////
  //  Calculate individual likelihood  //
  ///////////////////////////////////////

  for(i in 1:N){

    /////////////////////////////////////////////////////////////////
    // pSurvEntry: calculate probability of survival to entry      //
    //             into the cohort (account for left truncation).  //
    /////////////////////////////////////////////////////////////////

    if(lefttrunc_expose_DUR[i] > 0){
      log_psurventry <- log_sum_exp(diagonal(block(log_cumavoid_time_age, lefttrunc_expose_tIDX[i], lefttrunc_expose_aIDX[i], lefttrunc_expose_DUR[i], lefttrunc_expose_DUR[i])) +
      	      	                    diagonal(block(log_incrate_time_age, lefttrunc_expose_tIDX[i], lefttrunc_expose_aIDX[i], lefttrunc_expose_DUR[i], lefttrunc_expose_DUR[i])) +
                                    calc_log_phivsurv(lefttrunc_expose_tIDX[i], lefttrunc_expose_aIDX[i], lefttrunc_expose_taoIDX[i], lefttrunc_expose_DUR[i], entry_tIDX[i],
				                      log_hivsurv_REVdur_a0, log_hivmxMID_dur_a0, log_artrr_MID, artstart_tIDX, logdt)) +
			logdt + log_natsurv_time_age[entry_tIDX[i], entry_aIDX[i]];
       log_psurventry <- log_sum_exp(log_psurventry, log_natsurv_time_age[entry_tIDX[i], entry_aIDX[i]] + log_cumavoid_time_age[entry_tIDX[i], entry_aIDX[i]]);
    } else
      log_psurventry <- 0.0;


    ///////////////////////////////////////////////////////////////////////
    // phivn: probability of surviving and remaining HIV- [entry, exit]  //
    ///////////////////////////////////////////////////////////////////////

    if(hivpos[i])
      log_phivn <- negative_infinity();
    else
     log_phivn <- log_natsurv_time_age[exit_tIDX[i], exit_aIDX[i]] + log_cumavoid_time_age[exit_tIDX[i], exit_aIDX[i]] + death[i]*log_natmx_time_age[exit_tIDX[i], exit_aIDX[i]];

    ////////////////////////////////////////////////////////////////////////
    // phivp: probability of converting HIV+ and survival status at exit  //
    ////////////////////////////////////////////////////////////////////////

    if(expose_DUR[i] > 0){
      log_phivp <- log_sum_exp(diagonal(block(log_cumavoid_time_age, exposestart_tIDX[i], exposestart_aIDX[i], expose_DUR[i], expose_DUR[i])) + 
    	                    diagonal(block(log_incrate_time_age, exposestart_tIDX[i], exposestart_aIDX[i], expose_DUR[i], expose_DUR[i])) +
			    calc_log_phivsurv(exposestart_tIDX[i], exposestart_aIDX[i], exposestart_taoIDX[i], expose_DUR[i], exit_tIDX[i],
			                      log_hivsurv_REVdur_a0, log_hivmxMID_dur_a0, log_artrr_MID, artstart_tIDX, logdt) +
			    death[i]*log(exp(calc_log_hivmx(exposestart_tIDX[i], exposestart_aIDX[i], exposestart_taoIDX[i], expose_DUR[i], exit_tIDX[i],
			                                    log_hivmx_REVdur_a0, log_artrr, artstart_tIDX)) +
					 exp(log_natmx_time_age[exit_tIDX[i], exit_aIDX[i]]))) +
                 logdt + log_natsurv_time_age[exit_tIDX[i], exit_aIDX[i]];

    } else
      log_phivp <- negative_infinity();

    increment_log_prob(log_sum_exp(log_phivp, log_phivn) - log_psurventry);
  }

  print(get_lp(), " ", coef_art);
}