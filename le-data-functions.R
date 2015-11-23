## Load dataset
load("../alpha-prepare-data/ALPHA-data.RData")
hiv <- subset(hiv, type %in% c("survey") & !is.na(result) & !is.na(testdate))


## Declare functions

prepare.interval.data <- function(sites, sexes, min.age, max.age, min.time, max.time, hivonly=FALSE, hivelig=FALSE, nohiv=FALSE){

  ## ################## ##
  ##  Prepare the data  ##
  ## ################## ##

  ind <- subset(ind, site %in% sites & sex %in% sexes)
  res <- subset(res, site %in% sites & sex %in% sexes)
  if(nohiv)
    hiv <- subset(hiv, FALSE)
  else
    hiv <- subset(hiv, site %in% sites & sex %in% sexes)

  dat <- ind

  ## ################### ##
  ##  Truncate HIV data  ##
  ## ################### ##
  
  ## truncate HIV data range: only use HIV tests that occurred during age and time limits
  hiv$testage <- hiv$testdate - hiv$dob
  hiv <- subset(hiv, testdate >= min.time & testdate <= max.time & testage >= min.age & testage <= max.age)

  ## ################################ ##
  ##  Calculate observation episodes  ##
  ## ################################ ##
  
  if(!hivonly){
    ## identify observation intervals
    entry <- res[order(res$entry), c("id", "entry")]
    entry <- entry[!duplicated(entry$id),]
    exit <- res[order(res$exit, res$death, decreasing=TRUE), c("id", "exit", "death")]  # if 2 episodes exist ending on the same day, choose one with a death
    exit <- exit[!duplicated(exit$id),]
    
    ## create individual observation episodes
    dat <- merge(dat, entry)
    dat <- merge(dat, exit)

    if(hivelig){
      ## truncate to firsttest if hivelig=TRUE
      firsttest <- with(hiv, tapply(testdate, id, min))
      firsttest <- data.frame(id=as.integer(names(firsttest)), firsttest=firsttest)
      dat <- merge(dat, firsttest) # keep only people with HIV test
      dat$entry <- pmax(dat$entry, dat$firsttest)
      dat$firsttest <- NULL
    }
    
    ## truncate observation episodes
    dat <- subset(dat, entry <= max.time & exit >= min.time)  # eliminate episodes that do not intersect [min.time, max.time]
    dat$entry[dat$entry < min.time] <- min.time
    dat$death[dat$exit > max.time] <- 0
    dat$exit[dat$exit > max.time] <- max.time
    
    dat <- subset(dat, (entry-dob) <= max.age & (exit-dob) >= min.age)
    entryage <- dat$entry - dat$dob
    exitage <- dat$exit - dat$dob
    
    dat$death[exitage > max.age] <- 0
    entryage[entryage < min.age] <- min.age
    exitage[exitage > max.age] <- max.age
    
    dat$entry <- dat$dob + entryage
    dat$exit <- dat$dob + exitage
  } else { # hivonly = TRUE
    entry <- with(hiv, tapply(testdate, id, min))
    exit <- with(hiv, tapply(testdate, id, max))
    entry <- data.frame(id=as.integer(names(entry)), entry=entry)
    exit <- data.frame(id=as.integer(names(exit)), exit=exit)

    dat <- merge(dat, entry)
    dat <- merge(dat, exit)
    dat$death <- 0
  }
    
    
  ## ############## ##
  ##  Add HIV data  ##
  ## ############## ##

  ## create individual HIV test intervals
  lastneg <- with(subset(hiv, result == FALSE), tapply(testdate, id, max))
  firstpos <- with(subset(hiv, result == TRUE), tapply(testdate, id, min))
  lastneg <- data.frame(id=as.integer(names(lastneg)), lastneg=lastneg)
  firstpos <- data.frame(id=as.integer(names(firstpos)), firstpos=firstpos)

  ## merge HIV data
  dat <- merge(dat, lastneg, all.x=TRUE)
  dat <- merge(dat, firstpos, all.x=TRUE)

  dat$firstpos[dat$firstpos > dat$exit] <- NA
  dat$lastneg[dat$lastneg < dat$entry] <- NA
  dat$firstpos <- pmax(dat$firstpos, dat$entry)
  dat$lastneg <- pmin(dat$lastneg, dat$exit)
  
  attr(dat, "min.time") <- min.time
  attr(dat, "max.time") <- max.time
  attr(dat, "min.age") <- min.age
  attr(dat, "max.age") <- max.age

  return(dat)
}

discretise.cohort.data <- function(dat, dt){

  ## round dates to time steps [NOTE: think more about this -- should I take the floor instead?]
  dat$dobTS <- round(dat$dob / dt)
  dat$entryTS <- round(dat$entry / dt)
  dat$exitTS <- round(dat$exit / dt)
  dat$lastnegTS <- round(dat$lastneg / dt)
  dat$firstposTS <- round(dat$firstpos / dt)

  min.timeTS <- round(attr(dat, "min.time") / dt)
  max.timeTS <- round(attr(dat, "max.time") / dt)
  min.ageTS <- round(attr(dat, "min.age") / dt)
  max.ageTS <- round(attr(dat, "max.age") / dt)

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
  firstposageTS <- dat$firstposTS - dat$dobTS
  dat <- subset(dat, is.na(firstposageTS) | firstposageTS > min.ageTS)  # exclude people positive at first observation


  ## calculate incidence exposure intervals
  ## [exposestart, exposeend] defines the interval over which the person could have seroconverted
  dat$exposestartTS <- pmax(dat$lastnegTS, min.timeTS, dat$dobTS + min.ageTS, na.rm=TRUE)  ## ASSUMPTION: incidence doesn't occur before min.time or min.age
  dat$exposeendTS <- pmin(dat$firstposTS, dat$exitTS, na.rm=TRUE) - 1L  # -1 for must convert before observed positive
  dat$hivpos <- as.integer(!is.na(dat$firstposTS))  # flag for whether person is HIV+ at end of exposure period

  zeroexposehivp <- dat$hivpos & dat$exposestartTS > dat$exposeendTS  # if hivpos and start & end round to same time, assume infection occurred in previous time step. 
  dat$exposestartTS[zeroexposehivp] <- dat$exposestartTS[zeroexposehivp] - 1L 
  ## TODO: reconsider this later -- should date rounding be done differently?

  ## calculate array indices
  dat$entry.aIDX <- as.integer(dat$entryTS - dat$dobTS - min.ageTS) + 1L
  dat$exit.aIDX <- as.integer(dat$exitTS - dat$dobTS - min.ageTS) + 1L
  dat$exposestart.aIDX <- as.integer(dat$exposestartTS - dat$dobTS - min.ageTS) + 1L
  dat$exposeend.aIDX <- as.integer(dat$exposeendTS - dat$dobTS - min.ageTS) + 1L

  dat$entry.tIDX <- as.integer(dat$entryTS - min.timeTS) + 1L
  dat$exit.tIDX <- as.integer(dat$exitTS - min.timeTS) + 1L
  dat$exposestart.tIDX <- as.integer(dat$exposestartTS - min.timeTS) + 1L
  dat$exposeend.tIDX <- as.integer(dat$exposeendTS - min.timeTS) + 1L

  dat$expose.DUR <- dat$exposeend.tIDX - dat$exposestart.tIDX + 1L  # length(start:end) = end-start+1

  return(dat);  
}

aggregate.cohort.data <- function(dat){

  ## ############################### ##
  ##  Create aggregated cohort data  ##
  ## ############################### ##

  dat$cIDX <- dat$entry.tIDX - dat$entry.aIDX

  aggr <- dat[c("cIDX", "entry.tIDX", "exit.tIDX", "exposestart.tIDX", "exposeend.tIDX", "hivpos", "death")]
  aggr$nrepl <- 1

  ## left truncation intervals (psurventry)
  ltaggr <- data.frame(cIDX             = aggr$cIDX,
                       entry.tIDX       = pmax(aggr$cIDX+1, 1),
                       exit.tIDX        = aggr$entry.tIDX,
                       exposestart.tIDX = pmax(aggr$cIDX+1, 1),
                       exposeend.tIDX   = aggr$entry.tIDX-1,
                       hivpos           = 0,
                       death            = 0,
                       nrepl            = -1)  # subtract from each ll
  ltaggr <- subset(ltaggr, exposestart.tIDX <= exposeend.tIDX) # contributes 0 to likelihood if no exposure before entry (first observed at age.min)
  
  aggr <- rbind(aggr, ltaggr)
  aggr$entry.tIDX <- NULL

  aggr <- aggregate(nrepl~cIDX+exit.tIDX+exposestart.tIDX+exposeend.tIDX+hivpos+death, sum, data=aggr)
  aggr <- subset(aggr, nrepl != 0)
  aggr <- aggr[order(aggr$cIDX, aggr$exit.tIDX, -aggr$exposeend.tIDX),]
    
  ## for each cohort: minimum exposure time, max exposure time, number exit times, number individuals
  cohdat <- data.frame(t(sapply(split(aggr, aggr$cIDX),
                                function(coh) with(coh, c(coh_cIDX           = cIDX[1],
                                                          coh_minexpose_tIDX = min(exposestart.tIDX),
                                                          coh_maxexpose_tIDX = max(exposeend.tIDX),
                                                          coh_nexit          = length(unique(exit.tIDX)),
                                                          coh_ndat           = nrow(coh))))))

  ## for each (cohIDX, exitIDX): min exposure time for any future exit, max exposure time for any future exit
  exitdat <- lapply(split(aggr, aggr$cIDX), function(coh)
    data.frame(t(sapply(unique(coh$exit.tIDX), function(tidx) with(subset(coh, exit.tIDX >= tidx),
                                                                   c(exdat_cIDX           = cIDX[1],
                                                                     exdat_tIDX           = tidx,
                                                                     exdat_minexpose_tIDX = min(exposestart.tIDX),
                                                                     exdat_maxexpose_tIDX = min(tidx-1L, max(exposeend.tIDX)),
                                                                     exdat_ndat           = sum(exit.tIDX==tidx)))))))
  exitdat <- do.call(rbind, exitdat)

  return(list(cohdat=cohdat,
              exitdat=exitdat,
              aggr=aggr))
}
