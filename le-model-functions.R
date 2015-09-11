library(splines)
library(rstan)
library(expm)


#####################################
####  HIV survival distribution  ####
#####################################

hivsurv <- function(tao, a0, param){
  ## tao: time (years) since infection
  ## a0:  age (years) at infection

  ## Note: from Bellan et al. Lancet 2013, based on CASCADE Lancet 2000.
  shape <- 2.3
  scale <- 166/(shape*a0^0.53)
  return(exp(-(tao/scale)^shape))
}

hivsurvhaz <- function(tao, a0, param){
  shape <- 2.3
  scale <- 166/(shape*a0^0.53)
  return(shape/scale*(tao/scale)^(shape-1))
}


########################
####  Prepare data  ####
########################
                                   

prepare.stan.data <- function(sites = NULL, sexes = NULL, dat = NULL, dt = 0.1,
                              min.age = 15.0, max.age = 60.0,
                              min.time = 1980.5, max.time = 2011.5,
                              natmxstart.time, artstart.time,
                              ## nk.time = 7, nk.age = 7, nk.natmx = 5, nk.art = 5,
                              k.dt = 5, nk.art=5,
                              pen.ord.incrate=1L, pen.ord.natmx.time=1L, pen.ord.natmx.age=1L, pen.ord.art=1L,
                              nsamp=NULL, hivonly=FALSE, hivelig=FALSE){

  ## hivonly: if TRUE, only use HIV test data, don't use any residency episode data (no mortality).
  ## hivelig: indicates individuals only included if they have some HIV status information, so left truncate
  ##          at first HIV status information (inclusion conditional on survival to that point).

  if(is.null(dat)){
    dat <- prepare.interval.data(sites, sexes, min.age, max.age, min.time, max.time, hivonly, hivelig)
  }

  ## Select sub-sample 
  if(is.null(nsamp))
    samp <- 1:nrow(dat)  # all data
  else if (length(nsamp) == 1)
    samp <- sample(nrow(dat), nsamp)
  else
    samp <- nsamp  # specific indices specified

  dat <- dat[samp,]


  ## Discretise the dataset
  dat <- discretise.cohort.data(dat, dt=0.2)
  
  ##  Create aggregated cohort data for likelihood
  aggrdat <- aggregate.cohort.data(dat)
  aggr <- aggrdat$aggr
  exitdat <- aggrdat$exitdat
  cohdat <- aggrdat$cohdat



  min.timeTS <- round(min.time / dt)
  max.timeTS <- round(max.time / dt)
  min.ageTS <- round(min.age / dt)
  max.ageTS <- round(max.age / dt)


  artstart.timeTS <- round(artstart.time / dt)
  natmxstart.timeTS <- round(natmxstart.time / dt)
  
  artstart.tIDX <- as.integer(artstart.timeTS - min.timeTS) + 1L
  natmxstart.tIDX <- as.integer(natmxstart.timeTS - min.timeTS) + 1L


  ## ###################### ##
  ##  Prepare spline model  ##
  ## ###################### ##

  ## Model based on equally spaced knots through domain at interval k.dt

  x.time <- min.timeTS:max.timeTS*dt
  x.age <- min.ageTS:max.ageTS*dt
  STEPS_time=length(x.time)
  STEPS_age=length(x.age)


  ## incidence model
  k.incrate.time <- k.dt*(floor(min.time / k.dt) - 3L):(ceiling(max.time / k.dt) + 3L)
  k.incrate.age <- k.dt*(floor(min.age / k.dt) - 3L):(ceiling(max.age / k.dt) + 3L)

  nk_incrate_time <- length(k.incrate.time)-4L
  nk_incrate_age <- length(k.incrate.age)-4L

  X_incrate_time <- splineDesign(k.incrate.time, x.time, outer.ok=TRUE)
  Xmid_incrate_time <- splineDesign(k.incrate.time, x.time[-1]-dt/2, outer.ok=TRUE)
  X_incrate_age <- splineDesign(k.incrate.age, x.age, outer.ok=TRUE)
  Xmid_incrate_age <- splineDesign(k.incrate.age, x.age[-1]-dt/2)

  D_incrate_time <- diff(diag(nk_incrate_time), diff=pen.ord.incrate)
  D_incrate_age <- diff(diag(nk_incrate_age), diff=pen.ord.incrate)

  Pcar_prec_incrate <- matrix(0, nk_incrate_time*nk_incrate_age, nk_incrate_time*nk_incrate_age)
  diag(Pcar_prec_incrate[-1,]) <- rep(rep(c(0, -1), c(1, nk_incrate_time-1)), nk_incrate_age)[-1]
  diag(Pcar_prec_incrate[-(1:nk_incrate_time),]) <- -1
  Pcar_prec_incrate <- Pcar_prec_incrate+t(Pcar_prec_incrate)
  diag(Pcar_prec_incrate) <- -rowSums(Pcar_prec_incrate)

  Pcar_prec_incrate <- Pcar_prec_incrate %^% pen.ord.incrate


  ## non-HIV mortality model
  x.natmx <- natmxstart.timeTS:max.timeTS*dt

  k.natmx.time <- k.dt*(floor(natmxstart.time / k.dt) - 3L):(ceiling(max.time / k.dt) + 3L)
  k.natmx.age <- k.incrate.age

  nk_natmx_time <- length(k.natmx.time)-4L
  nk_natmx_age <- length(k.natmx.age)-4L

  X_natmx_time <- splineDesign(k.natmx.time, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx), outer.ok=TRUE)
  Xmid_natmx_time <- splineDesign(k.natmx.time, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx[-1]-dt/2))
  X_natmx_age <- splineDesign(k.natmx.age, x.age, outer.ok=TRUE)
  Xmid_natmx_age <- splineDesign(k.natmx.age, x.age[-1]-dt/2)

  D_natmx_time<- diff(diag(nk_natmx_time), diff=pen.ord.natmx.time)
  D_natmx_age <- diff(diag(nk_natmx_age), diff=pen.ord.natmx.age)

  
  ## ART model

  ## x.art <- artstart.timeTS:max.timeTS*dt
  ## X.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splineDesign(k.art, x.art, outer.ok=TRUE))
  ## Xmid.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splineDesign(k.art, x.art[-1]-dt/2))
  ## P.art <- diff(diag(nk.art), diff=1)
  
  D_art <- diff(diag(STEPS_time - artstart.tIDX), diff=pen.ord.art)

  ## ##################################### ##
  ##  Calculate HIV survival lookup table  ##
  ## ##################################### ##

  log_hivmx_dur_a0 <- log(outer(1:(max.timeTS-min.timeTS)*dt - dt/2, min.ageTS:(max.ageTS-1L)*dt + dt/2, hivsurvhaz, param=NULL))
  log_hivsurv_dur_a0 <- log(outer(1:(max.timeTS-min.timeTS)*dt - dt/2, min.ageTS:(max.ageTS-1L)*dt + dt/2, hivsurv, param=NULL))
  log_hivmxMID_dur_a0 <- log(-diff(rbind(0, log_hivsurv_dur_a0))/dt)


  hivmx_dur_a0 <- exp(log_hivmx_dur_a0)
  hivsurv_dur_a0 <- exp(log_hivsurv_dur_a0)
  hivmxMID_dur_a0 <- exp(log_hivmxMID_dur_a0)

    
  ## ######################## ##
  ##  Create Stan input data  ##
  ## ######################## ##

  stan.data <- list(dt                    = dt,
                    STEPS_time            = STEPS_time,
                    STEPS_age             = STEPS_age,
                    artstart_tIDX         = artstart.tIDX,
                    ## INDIVIDUAL DATA  ##
                    N                     = nrow(dat),
                    id                    = dat$id,
                    entry_tIDX            = dat$entry.tIDX,
                    exit_tIDX             = dat$exit.tIDX,
                    exposestart_tIDX      = dat$exposestart.tIDX,
                    entry_aIDX            = dat$entry.aIDX,
                    exit_aIDX             = dat$exit.aIDX,
                    exposestart_aIDX      = dat$exposestart.aIDX,
                    expose_DUR            = dat$expose.DUR,
                    hivpos                = dat$hivpos,
                    death                 = dat$death,
                    ## COHORT DATA ##
                    NCOH                  = nrow(cohdat),
                    coh_cIDX              = cohdat$coh_cIDX,
                    coh_minexpose_tIDX    = cohdat$coh_minexpose_tIDX,
                    coh_maxexpose_tIDX    = cohdat$coh_maxexpose_tIDX,
                    coh_nexit             = cohdat$coh_nexit,
                    coh_ndat              = cohdat$coh_ndat,
                    ## EXIT DATA ##
                    NEXIT                 = nrow(exitdat),
                    exdat_cIDX            = exitdat$exdat_cIDX,
                    exdat_tIDX            = exitdat$exdat_tIDX,
                    exdat_minexpose_tIDX  = exitdat$exdat_minexpose_tIDX,
                    exdat_maxexpose_tIDX  = exitdat$exdat_maxexpose_tIDX,
                    exdat_ndat            = exitdat$exdat_ndat,
                    ## AGGR DATA ##
                    NAGGR                 = nrow(aggr),
                    aggr_cIDX             = aggr$cIDX,
                    aggr_exit_tIDX        = aggr$exit.tIDX,
                    aggr_exposestart_tIDX = aggr$exposestart.tIDX,
                    aggr_exposeend_tIDX   = aggr$exposeend.tIDX,
                    aggr_death            = aggr$death,
                    aggr_hivpos           = aggr$hivpos,
                    aggr_nrepl            = aggr$nrepl,
                    ## MODEL PARAMETERS ##
                    x_time                = x.time,
                    x_age                 = x.age,
                    ## x_art                 = x.art,
                    x_natmx               = x.natmx,
                    ## incidence model
                    nk_incrate_time       = nk_incrate_time,
                    nk_incrate_age        = nk_incrate_age,
                    k_incrate_time        = k.incrate.time,
                    k_incrate_age         = k.incrate.age,
                    X_incrate_time        = X_incrate_time,
                    Xmid_incrate_time     = Xmid_incrate_time,
                    X_incrate_age         = X_incrate_age,
                    Xmid_incrate_age      = Xmid_incrate_age,
                    D_incrate_time        = D_incrate_time,
                    D_incrate_age         = D_incrate_age,
                    pen_ord_incrate       = pen.ord.incrate,
                    Pcar_prec_incrate     = Pcar_prec_incrate,
                    ## non-HIV mortality model
                    nk_natmx_time         = nk_natmx_time,
                    nk_natmx_age          = nk_natmx_age,
                    k_natmx_time          = k.natmx.time,
                    k_natmx_age           = k.natmx.age,
                    X_natmx_time          = X_natmx_time,
                    Xmid_natmx_time       = Xmid_natmx_time,
                    X_natmx_age           = X_natmx_age,
                    Xmid_natmx_age        = Xmid_natmx_age,
                    pen_ord_natmx_time    = pen.ord.natmx.time,
                    pen_ord_natmx_age     = pen.ord.natmx.age,
                    D_natmx_time          = D_natmx_time,
                    D_natmx_age           = D_natmx_age,
                    fixcoef_natmx_time    = as.integer(nk_natmx_time/2),
                    fixcoef_natmx_age     = as.integer(nk_natmx_age/2),
                    ## ART model
                    ## nk_art                = nk.art,
                    ## k_art                 = k.art,
                    ## X_art                 = X.art,
                    ## Xmid_art              = Xmid.art,
                    pen_ord_art           = pen.ord.art,
                    D_art                 = D_art,
                    ## HIV survival model
                    log_hivmx_dur_a0      = log_hivmx_dur_a0,
                    log_hivsurv_dur_a0    = log_hivsurv_dur_a0,
                    log_hivmxMID_dur_a0   = log_hivmxMID_dur_a0,
                    hivmx_dur_a0          = hivmx_dur_a0,
                    hivsurv_dur_a0        = hivsurv_dur_a0,
                    hivmxMID_dur_a0       = hivmxMID_dur_a0)
}



##################################
####                          ####
####    Analysis functions    ####
####                          ####
##################################

invlogit <- function(x) exp(x)/(1+exp(x))

convert.stan.params <- function(par, stand){

  nk_time <- stand$nk_time
  nk_age <- stand$nk_age
  nk_natmx <- stand$nk_natmx
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

  incrateMID_time_age <- exp(stand$Xmid_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_age))
  cumavoid_time_age <- exp(-stand$dt*diagCumSum(incrateMID_time_age))
  cumavoidMID_time_age <- cumavoid_time_age[1:(stand$STEPS_time-1), 1:(stand$STEPS_age-1)]*exp(-stand$dt/2*incrateMID_time_age)

  natmx_time_age <- exp(stand$X_natmx %*% param$coef_natmx_time_age %*% t(stand$X_age))
  natsurv_time_age <- exp(-stand$dt*diagCumSum(exp(stand$Xmid_natmx %*% param$coef_natmx_time_age %*% t(stand$Xmid_age))))

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
  log_incrate_time_age <- stand$X_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_age)

  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  setNames(cumincid.period, stand$x_time)
  return(cumincid.period)
}

cumnatmort.period <- function(param, stand){
  log_natmx_time_age <- stand$X_natmx %*% param$coef_natmx_time_age %*% t(stand$Xmid_age)

  cumnatmort.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_natmx_time_age)))
  setNames(cumnatmort.period, stand$x_time)
  return(cumnatmort.period)
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



#### Calculate likelihood

Rcalc_ll_coh <- function(stand, modpred){
  sum(calc_ll_coh(stand$entry_tIDX, stand$entry_aIDX, stand$exit_tIDX, stand$exit_aIDX,
                  stand$exposestart_tIDX, stand$exposestart_aIDX, stand$expose_DUR, stand$death, stand$hivpos,
                  modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                  stand$hivsurv_dur_a0, stand$hivmx_dur_a0, stand$hivmxMID_dur_a0,
                  modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                  modpred$natsurv_time_age, modpred$natmx_time_age,
                  stand$dt))
}

Rcalc_ll_coh_aggr <- function(stand, modpred){
  calc_ll_coh_aggr(stand$aggr_cIDX, stand$aggr_exit_tIDX, stand$aggr_exposestart_tIDX, stand$aggr_exposeend_tIDX, stand$aggr_death, stand$aggr_hivpos, stand$aggr_nrepl,
                   modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                   stand$hivsurv_dur_a0, stand$hivmx_dur_a0, stand$hivmxMID_dur_a0,
                   modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                   modpred$natsurv_time_age, modpred$natmx_time_age,
                   stand$dt)
}

Rcalc_ll_cohexit <- function(stand, modpred){
  calc_ll_cohexit(stand$coh_cIDX, stand$coh_minexpose_tIDX, stand$coh_maxexpose_tIDX, stand$coh_nexit,
                  stand$exdat_tIDX, stand$exdat_minexpose_tIDX, stand$exdat_maxexpose_tIDX, stand$exdat_ndat,
                  stand$aggr_exposestart_tIDX, stand$aggr_exposeend_tIDX, stand$aggr_death, stand$aggr_hivpos, stand$aggr_nrepl,
                  modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                  stand$hivsurv_dur_a0, stand$hivmx_dur_a0, stand$hivmxMID_dur_a0,
                  modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                  modpred$natsurv_time_age, modpred$natmx_time_age,
                  stand$dt)
}
