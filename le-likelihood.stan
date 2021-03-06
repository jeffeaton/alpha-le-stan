  ///////////////////////////////////////////////////////
  //  Construct incidence rate and mortality matrices  //
  ///////////////////////////////////////////////////////

  {
    matrix[STEPS_time-1, STEPS_age-1] incrateMID_time_age;
    matrix[STEPS_time, STEPS_age] cumavoid_time_age;
    matrix[STEPS_time-1, STEPS_age-1] cumavoidMID_time_age;
    matrix[STEPS_time, STEPS_age] natmx_time_age;
    matrix[STEPS_time, STEPS_age] natsurv_time_age;
    vector[STEPS_time] artrr;
    vector[STEPS_time-1] artrr_MID;
    
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
}
