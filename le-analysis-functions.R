##################################
####                          ####
####    Analysis functions    ####
####                          ####
##################################

invlogit <- function(x) exp(x)/(1+exp(x))

convert.stan.params <- function(par, stand){

  nk_time <- stand$nk_time
  nk_age <- stand$nk_age
  nk_natmx_time <- stand$nk_natmx_time
  nk_art <- stand$nk_art
  
  coef_incrate_time_age.idx <- 1:(nk_time*nk_age)
  coef_natmx_time.idx <- nk_time*nk_age + 1:nk_natmx
  coef_natmx_age.idx <- nk_time*nk_age + nk_natmx + 1:(nk_age-1L)
  coef_art.idx <- nk_time*nk_age + nk_natmx + (nk_age-1L) + 1:nk_art

  coef_incrate_time_age <- matrix(par[coef_incrate_time_age.idx], nk_time, nk_age)

  coef_natmx_time <- par[coef_natmx_time.idx]
  
  coef_natmx_age <- rep(0, nk_age)
  coef_natmx_age[-stand$fixcoef_age_idx] <- par[coef_natmx_age.idx]
    
  coef_natmx_time_age <- outer(coef_natmx_time, coef_natmx_age, "+")
  
  coef_art <- par[coef_art.idx]

  return(list(coef_incrate_time_age = coef_incrate_time_age,
              coef_natmx_time_age = coef_natmx_time_age,
              coef_natmx_time = coef_natmx_time,
              coef_natmx_age = coef_natmx_age,
              coef_art = coef_art))
}

create.param.list <- function(stanfit){
  param <- extract(stanfit)
  param <- lapply(seq_along(param$lp__), function(ii) list(coef_incrate_time_age   = param$coef_incrate_time_age[ii,,],
                                                           coef_natmx_time         = param$coef_natmx_time[ii,],
                                                           coef_natmx_age          = param$coef_natmx_age[ii,],
                                                           coef_natmx_time_age     = outer(param$coef_natmx_time[ii,], param$coef_natmx_age[ii,], "+"),
                                                           ## coef_art                = param$coef_art[ii,],
                                                           dt_log_artrr            = param$dt_log_artrr[ii,],
                                                           sigma2_incrate_time_age = param$sigma2_incrate_time_age[ii],
                                                           sigma2_natmx_time       = param$sigma2_natmx_time[ii],
                                                           sigma2_natmx_age        = param$sigma2_natmx_age[ii],
                                                           sigma2_art              = param$sigma2_art[ii]))
  return(param)
}

create.modpred <- function(param, stand){

  incrateMID_time_age <- exp(stand$Xmid_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age))
  cumavoid_time_age <- exp(-stand$dt*diagCumSum(incrateMID_time_age))
  cumavoidMID_time_age <- cumavoid_time_age[1:(stand$STEPS_time-1), 1:(stand$STEPS_age-1)]*exp(-stand$dt/2*incrateMID_time_age)

  natmx_time_age <- exp(stand$X_natmx_time %*% param$coef_natmx_time_age %*% t(stand$X_natmx_age))
  natsurv_time_age <- exp(-stand$dt*diagCumSum(exp(stand$Xmid_natmx_time %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age))))

  ## artrr <- invlogit(stand$X_art %*% param$coef_art)
  ## artrr_MID <- invlogit(stand$Xmid_art %*% param$coef_art)
  ## artrr[1:(stand$artstart_tIDX-1)] <- 1.0
  ## artrr_MID[1:(stand$artstart_tIDX-1)] <- 1.0

  log_artrr <- stand$dt*cumsum(param$dt_log_artrr);
  artrr <- rep(1.0, stand$STEPS_time)
  artrr_MID <- rep(1.0, stand$STEPS_time-1L)
  artrr[(stand$artstart_tIDX+1L):stand$STEPS_time] <- exp(log_artrr)
  artrr_MID[stand$artstart_tIDX:(stand$STEPS_time-1L)] <- exp(log_artrr - stand$dt/2*param$dt_log_artrr)

  return(list(incrateMID_time_age = incrateMID_time_age,
              cumavoid_time_age = cumavoid_time_age,
              cumavoidMID_time_age = cumavoidMID_time_age,
              natmx_time_age = natmx_time_age,
              natsurv_time_age = natsurv_time_age,
              artrr = artrr,
              artrr_MID = artrr_MID))
}

cumincid.period <- function(param, stand){
  log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age)

  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  return(setNames(cumincid.period, stand$x_time))
}

cumnatmort.period <- function(param, stand){
  log_natmx_time_age <- stand$X_natmx_time %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age)

  cumnatmort.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_natmx_time_age)))
  return(setNames(cumnatmort.period, stand$x_time))
}

prev <- function(tidx, aidx, modpred, stand){
  exposeDUR <- min(tidx, aidx)-1L
  phivp <- calc_phivp(tidx, aidx, tidx-exposeDUR, aidx-exposeDUR, exposeDUR, 0,
                      modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                      stand$hivsurv_dur_a0, stand$hivmx_dur_a0, stand$hivmxMID_dur_a0, modpred$artrr, modpred$artrr_MID,
                      stand$artstart_tIDX, modpred$natsurv_time_age, modpred$natmx_time_age, stand$dt)
  phivn <- calc_phivn(tidx, aidx, 0, 0, modpred$cumavoid_time_age, modpred$natsurv_time_age, modpred$natmx_time_age)
  return(phivp/(phivp+phivn))
}


Rcreate_phivp_mat <- function(modpred, stand, art=TRUE){
  if(!art)
    stand$artstart_tIDX <- stand$STEPS_time
  create_phivp_mat(modpred$cumavoidMID_time_age, modpred$incrateMID_time_age, stand$hivsurv_dur_a0, stand$hivmxMID_dur_a0, modpred$artrr_MID, stand$artstart_tIDX, modpred$natsurv_time_age, stand$dt)
}
Rcreate_phivn_mat <- function(modpred){ return(modpred$cumavoid_time_age * modpred$natsurv_time_age)}

calc.psurv <- function(param, stand){
  modpred <- create.modpred(param, stand)
  phivn <- Rcreate_phivn_mat(modpred)
  phivp <- Rcreate_phivp_mat(modpred, stand)
  phivp.noart <- Rcreate_phivp_mat(modpred, stand, art=FALSE)
  psurv <- phivn+phivp
  psurv.noart <- phivn+phivp.noart
  psurv.nohiv <- modpred$natsurv_time_age
  list(psurv=psurv, psurv.noart=psurv.noart, psurv.nohiv=psurv.nohiv, prev=phivp/psurv, prev.noart=phivp.noart/psurv.noart)
}

calc.le <- function(psurv, dt){
  psurv.last <- psurv[-nrow(psurv), -ncol(psurv)]
  psurv.curr <- psurv[-1,-1]
  dt*colSums(apply(1.0 - (psurv.last-psurv.curr)/psurv.last, 1, cumprod))
}


###############################
####  Aggregate functions  ####
###############################

library(parallel)

calc.45q15 <- function(psurv, dt){
  aidx <- seq_len(45/dt) # assumes psurv starts at age 15
  psurv.last <- psurv[-nrow(psurv), -ncol(psurv)]
  psurv.curr <- psurv[-1,-1]
  qx <- (psurv.last[,aidx]-psurv.curr[,aidx])/psurv.last[,aidx]
  1.0 - exp(rowSums(log(1.0-qx)))
}

add.le <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(le       = sapply(lapply(psurvobj, "[[", "psurv"), calc.le, mod$stand$dt),
       le.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.le, mod$stand$dt),
       le.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.le, mod$stand$dt))
}

add.45q15 <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(q4515       = sapply(lapply(psurvobj, "[[", "psurv"), calc.45q15, mod$stand$dt),
       q4515.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.45q15, mod$stand$dt),
       q4515.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.45q15, mod$stand$dt))
}

add.mx <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(le       = sapply(lapply(psurvobj, "[[", "psurv"), calc.le, mod$stand$dt),
       le.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.le, mod$stand$dt),
       le.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.le, mod$stand$dt),
       q4515       = sapply(lapply(psurvobj, "[[", "psurv"), calc.45q15, mod$stand$dt),
       q4515.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.45q15, mod$stand$dt),
       q4515.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.45q15, mod$stand$dt))
}

