library(splines)
library(rstan)

#####################################
####  HIV survival distribution  ####
#####################################

hivsurv <- function(tao, a0, param){
  ## tao: time (years) since infection
  ## a0:  age (years) at infection

  ## Note: from Bellan et al. Lancet 2013, based on CASCADE Lancet 2000. 
  k <- 2.3
  Lambda <- 166/(k*a0^0.53)
  return(exp(-(tao/Lambda)^k))
}

hivsurvhaz <- function(tao, a0, param){
  k <- 2.3
  Lambda <- 166/(k*a0^0.53)
  return(k/Lambda*(tao/Lambda)^(k-1))
}


########################
####  Prepare data  ####
########################

prepare.data <- function(dat, dt = 0.1, min.age = 15.0, max.age = 85.0, min.time = 1980.5, max.time = 2011.5, natmxstart.time, artstart.time, nk.time = 7, nk.age = 7, nk.natmx = 5, nk.art = 5, nsamp=NULL){

  ## ####### ARGUMENTS ###########
##   dt <- 0.2
##   min.age <- 15.0
##   max.age <- 100
##   min.time <- 1980.5
## max.time <- 2013.0
## artstart.time <- 2005.5
## natmxstart.time <- 1990.0
##   nk.time <- 7
## nk.age <- 7
## nk.art <- 5
## nk.natmx <- 5
  ## ############################
  
  ## round dates to time steps [NOTE: think more about this -- should I take the floor instead?]
  dat$dobTS <- round(dat$dob / dt)
  dat$entryTS <- round(dat$entry / dt)
  dat$exitTS <- round(dat$exit / dt)
  ## dat$firsttestTS <- round(dat$firsttest / dt)
  dat$lastnegTS <- round(dat$lastneg / dt)
  dat$firstposTS <- round(dat$firstpos / dt)
  
  min.timeTS <- round(min.time / dt) ## ifelse(is.null(min.time), min(dat$firsttestTS), round(min.time / dt))
  max.timeTS <- round(max.time / dt) ## ifelse(is.null(max.time), max(dat$lastnegTS, dat$firstposTS, na.rm=TRUE), round(max.time / dt))
  artstart.timeTS <- round(artstart.time / dt)
  natmxstart.timeTS <- round(natmxstart.time / dt)
  min.ageTS <- round(min.age / dt)
  max.ageTS <- round(max.age / dt)
  
  ## dat <- subset(dat, entryTS < exitTS) # eliminate people with exitTS - entryTS = 0   !!!! REVISE THIS: if have HIV information or death, then includes info

  dat <- subset(dat, is.na(lastnegTS) | is.na(firstposTS) | lastnegTS <= firstposTS) # omit inconsistent HIV data
  
  
  ## age censoring   [NOTE: slight bias by censoring positive observations but not negative person-time]
  dat <- subset(dat, exitTS - dobTS >= min.ageTS)
  dat <- subset(dat, entryTS - dobTS <= max.ageTS)
  
  entryageTS <- dat$entryTS - dat$dobTS
  exitageTS <- dat$exitTS - dat$dobTS
  lastnegageTS <- dat$lastnegTS - dat$dobTS
  firstposageTS <- dat$firstposTS - dat$dobTS
  
  entryageTS[entryageTS < min.ageTS] <- min.ageTS
  dat$entryTS <- dat$dobTS + entryageTS
  
  dat$death[exitageTS > max.ageTS] <- 0
  exitageTS[exitageTS > max.ageTS] <- max.ageTS
  dat$exitTS <- dat$dobTS + exitageTS

  dat$lastnegTS[lastnegageTS < min.ageTS] <- NA
  lastnegageTS[lastnegageTS > max.ageTS] <- max.ageTS  # bias: should be censored to lastneg before censor date
  dat$lastnegTS <- dat$dobTS + lastnegageTS

  firstposageTS[firstposageTS < min.ageTS] <- min.ageTS
  firstposTS <- dat$dobTS + firstposageTS
  dat$firstposTS[firstposageTS > max.ageTS] <- NA

  ## TEMPORARY: eventually incorporate P(positive at age 15) into model
  dat <- subset(dat, is.na(firstposageTS) | firstposageTS > min.ageTS)

  rm(entryageTS, exitageTS, lastnegageTS, firstposageTS)
  

  ## time censoring
  dat <- subset(dat, exitTS >= min.timeTS)
  dat <- subset(dat, entryTS <= max.timeTS)
  
  dat$death[dat$exitTS > max.timeTS] <- 0
  dat$exitTS[dat$exitTS > max.timeTS] <- max.timeTS
  
  dat$entryTS[dat$entryTS < min.timeTS] <- min.timeTS
  
  dat$lastnegTS[dat$lastnegTS < min.timeTS] <- NA
  dat$lastnegTS[dat$lastnegTS > max.timeTS] <- max.timeTS  # bias: should be censored to lastneg before censor date
  
  dat$firstposTS[dat$firstposTS < min.timeTS] <- min.timeTS
  dat$firstposTS[dat$firstposTS > max.timeTS] <- NA
  
  
  ## calculate incidence exposure intervals
  ## [exposestart, exposeend] defines the interval over which the person could have seroconverted
  dat$exposestartTS <- pmax(dat$lastnegTS, min.timeTS, dat$dobTS + min.ageTS, na.rm=TRUE)  ## ASSUMPTION: incidence doesn't occur before min.time or min.age
  dat$exposeendTS <- pmin(dat$firstposTS, dat$exitTS, na.rm=TRUE)
  dat$hivpos <- as.integer(!is.na(dat$firstposTS))  # flag for whether person is HIV+ at end of exposure period
  
  
  ## calculate array indices
  
  dat$entry.aIDX <- as.integer(dat$entryTS - dat$dobTS - min.ageTS) + 1L
  dat$exit.aIDX <- as.integer(dat$exitTS - dat$dobTS - min.ageTS) + 1L
  dat$exposestart.aIDX <- as.integer(dat$exposestartTS - dat$dobTS - min.ageTS) + 1L
  dat$exposeend.aIDX <- as.integer(dat$exposeendTS - dat$dobTS - min.ageTS) + 1L
  
  dat$entry.tIDX <- as.integer(dat$entryTS - min.timeTS) + 1L
  dat$exit.tIDX <- as.integer(dat$exitTS - min.timeTS) + 1L
  dat$exposestart.tIDX <- as.integer(dat$exposestartTS - min.timeTS) + 1L
  dat$exposeend.tIDX <- as.integer(dat$exposeendTS - min.timeTS) + 1L

  dat$expose.DUR <- dat$exposeend.tIDX - dat$exposestart.tIDX
  
  artstart.tIDX <- as.integer(artstart.timeTS - min.timeTS) + 1L
  natmxstart.tIDX <- as.integer(natmxstart.timeTS - min.timeTS) + 1L
  
  
  ## ###################### ##
  ##  Prepare spline model  ##
  ## ###################### ##
  
  x.time <- min.timeTS:max.timeTS*dt
  x.age <- min.ageTS:max.ageTS*dt
  x.art <- artstart.timeTS:max.timeTS*dt
  x.natmx <- natmxstart.timeTS:max.timeTS*dt
  
  time.dur <- max.timeTS * dt - min.timeTS * dt
  k.time <- seq(min.timeTS*dt - 3*time.dur/(nk.time-3), max.timeTS*dt + 3*time.dur/(nk.time-3), time.dur/(nk.time-3))
  
  age.dur <- (max.ageTS - min.ageTS)*dt
  k.age <- seq(min.ageTS*dt - 3*age.dur/(nk.age-3), max.ageTS*dt + 3*age.dur/(nk.age-3), age.dur/(nk.age-3))
  
  art.dur <- (max.timeTS - artstart.timeTS)*dt
  k.art <- seq(artstart.timeTS*dt - 3*art.dur/(nk.art-3), max.timeTS*dt + 3*art.dur/(nk.art-3), art.dur/(nk.art-3))
  
  natmx.dur <- (max.timeTS - natmxstart.timeTS)*dt
  k.natmx <- seq(natmxstart.timeTS*dt - 3*natmx.dur/(nk.natmx-3), max.timeTS*dt + 3*natmx.dur/(nk.natmx-3), natmx.dur/(nk.natmx-3))
  
  Xmid.time <- splineDesign(k.time, x.time[-1]-dt/2, outer.ok=TRUE)
  Xmid.age <- splineDesign(k.age, x.age[-1]-dt/2)
  Xmid.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splineDesign(k.art, x.art[-1]-dt/2))
  Xmid.natmx <- splineDesign(k.natmx, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx[-1]-dt/2))
  
  
  X.time <- splineDesign(k.time, x.time)
  X.age <- splineDesign(k.age, x.age)
  X.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splineDesign(k.art, x.art))
  X.natmx <- splineDesign(k.natmx, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx))
  
  P.time <- diff(diag(nk.time), diff=1)
  P.age <- diff(diag(nk.age), diff=1)
  P.art <- diff(diag(nk.art), diff=1)
  P.natmx <- diff(diag(nk.natmx), diff=1)
  
  STEPS_time=length(x.time)
  STEPS_age=length(x.age)

  ## ##################################### ##
  ##  Calculate HIV survival lookup table  ##
  ## ##################################### ##

  log_hivmx_REVdur_a0 <- log(outer((max.timeTS-min.timeTS):0*dt, min.ageTS:max.ageTS*dt, hivsurvhaz, param=NULL))
  log_hivsurv_REVdur_a0 <- log(outer((max.timeTS-min.timeTS):0*dt, min.ageTS:max.ageTS*dt, hivsurv, param=NULL))
  log_hivmxMID_dur_a0 <- log(outer(1:(max.timeTS-min.timeTS)*dt-dt/2, min.ageTS:max.ageTS*dt, hivsurvhaz, param=NULL))
  
  ## ######################## ##
  ##  Create Stan input data  ##
  ## ######################## ##

  
  if(is.null(nsamp))
    samp <- 1:nrow(dat)
  else if (length(nsamp) == 1)
    samp <- sample(nrow(dat), nsamp)
  else
    samp <- nsamp

  stan.data <- list(dt                    = dt,
                    STEPS_time            = STEPS_time,
                    STEPS_age             = STEPS_age,
                     artstart_tIDX         = artstart.tIDX,
                     ##
                    N                     = length(samp),
                    idno                  = dat$idno[samp],
                     entry_tIDX            = dat$entry.tIDX[samp],
                     exit_tIDX             = dat$exit.tIDX[samp],
                     exposestart_tIDX      = dat$exposestart.tIDX[samp],
                     entry_aIDX            = dat$entry.aIDX[samp],
                     exit_aIDX             = dat$exit.aIDX[samp],
                     exposestart_aIDX      = dat$exposestart.aIDX[samp],
                     expose_DUR            = dat$expose.DUR[samp],
                     hivpos                = dat$hivpos[samp],
                     death                 = dat$death[samp],
                     ##
                     ## log_incrate_time_age  = log_incrate_time_age,
                     ## log_cumavoid_time_age = log_cumavoid_time_age,
                     ## log_natmx_time_age    = log_natmx_time_age,
                     ## log_natsurv_time_age  = log_natsurv_time_age,
                     log_hivmx_REVdur_a0      = log_hivmx_REVdur_a0,
                     log_hivsurv_REVdur_a0    = log_hivsurv_REVdur_a0,
                     log_hivmxMID_dur_a0      = log_hivmxMID_dur_a0,
                     ## log_artrr_REV         = log_artrr_REV)
                     ##
                     nk_time               = nk.time,
                     nk_age                = nk.age,
                     nk_art                = nk.art,
                     nk_natmx              = nk.natmx,
                     fixcoef_age_idx       = 3L,
                     fixcoef_time_idx      = 4L,
                     X_time                = X.time,
                     X_age                 = X.age,
                     X_art                 = X.art,
                     X_natmx               = X.natmx,
                     Xmid_time             = Xmid.time,
                     Xmid_age              = Xmid.age,
                     Xmid_art              = Xmid.art,
                     Xmid_natmx            = Xmid.natmx,                     
                     P_time                = P.time,
                     P_age                 = P.age,
                     P_art                 = P.art,
                     P_natmx               = P.natmx,
                     ##
                     sigma2_incrate_time   = 2.0,
                     sigma2_incrate_age    = 2.0,
                     sigma2_incrate_time_age = 0.3,
                     sigma2_natmx_time       = 0.1,
                    sigma2_art              = 20)
}
