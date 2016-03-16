library(curl)
library(rstan)

catfiles <- function(..., path="", url=TRUE, string=TRUE){
  dots <- list(...)
  if(url)
    lines <- do.call(c, lapply(lapply(paste(path, dots, sep=""), curl), readLines))
  else
    lines <- do.call(c, lapply(paste(path, dots, sep=""), readLines))

  if(string)
    return(paste(lines, collapse="\n"))
  else
    return(lines)
}

compile.chains <- function(prefix, path=""){
  files <- system(paste("ls ", path, prefix, "*.RData", sep=""), inter=TRUE)
  print(paste(prefix, ": ", length(files), " chains found", sep=""))
  if(length(files) == 0) return(NULL);
  lapply(files, load, environment())
  sflist2stanfit(lapply(grep(prefix, ls(), value=TRUE), get, envir=environment()))
}

expose_file <- function(stanfile, path="", url=TRUE){
  if(url)
    code <- readLines(curl(paste(path, stanfile, sep="")))
  else
    code <- readLines(paste(path, stanfile, sep=""))
  expose_stan_functions(stan_model(model_code=paste(c(code, "\nmodel{}"), collapse="\n")))
  return(invisible())
}


## ## Example:
## sourceurl <- "https://raw.github.com/jeffeaton/alpha-le-stan/master/"
## mod <- catfiles("le-functions.stan", "le-data.stan", "le-model.stan", path=sourceurl)



create.cluster.scripts <- function(label, site, sex, chains, iter=1000, min.time, max.time, natmxstart.time, artstart.time,
                                   min.age = 15.0, max.age = 100.0, dt=0.1, k.dt=5, nk.art=5,
                                   cohortstart.time=min.time, cohortend.time=max.time,
                                   pen.ord=1L, pen.ord.incrate=pen.ord, pen.ord.natmx.time=pen.ord, pen.ord.natmx.age=pen.ord,
                                   pen.ord.art=pen.ord, hivonly=FALSE, hivelig=FALSE, max.tree.depth=10,
                                   model='isotrop', sourceurl="https://raw.github.com/jeffeaton/alpha-le-stan/master/", clust.path=clustpath, file.path=filepath){

  ## Create submit-jobs.sh if it doesn't already exist
  if(!file.exists(paste(file.path, "submit-jobs.sh", sep=""))){
    fileConn <- file(paste(file.path, "submit-jobs.sh", sep=""))
    writeLines("#!/bin/bash", fileConn)
    close(fileConn)
  }

  for(i in chains){
    fileConn <- file(paste(file.path, label, "-fit", i, ".R", sep=""))
    writeLines(c(paste("setwd(paste(Sys.getenv('HOME'), '/", clust.path, "', sep=''))", sep=""),
                 paste("sourceurl <- '", sourceurl, "'", sep=""),
                 ## "library(devtools)",
                 "source(paste(sourceurl, 'utilities.R', sep=''))",
                 "source(paste(sourceurl, 'le-data-functions.R', sep=''))",
                 "source(paste(sourceurl, 'le-model-functions.R', sep=''))",
                 paste("mod <- catfiles('le-functions.stan', 'le-data.stan', 'le-model_", model, ".stan', 'le-likelihood.stan', path=sourceurl, url=FALSE)", sep=""),
                 paste(label, ".stand <- prepare.stan.data(c('", site, "'), '", sex, "'",
                       ", min.time=", min.time, ", max.time=", max.time,
                       ", min.age=", min.age, ", max.age=", max.age,
                       ", natmxstart.time=", natmxstart.time, ", artstart.time=", artstart.time,
                       ", cohortstart.time=", cohortstart.time, ", cohortend.time=", cohortend.time,
                       ", pen.ord.incrate=", pen.ord.incrate, ", pen.ord.natmx.time=", pen.ord.natmx.time,
                       ", pen.ord.natmx.age=", pen.ord.natmx.age, ", pen.ord.art=", pen.ord.art, 
                       ", dt=", dt, ", k.dt=", k.dt, ", nk.art=", nk.art, ", hivonly=", hivonly, ", hivelig=", hivelig, ")",
                       sep=""),
                 paste(label, ".fit", i, " <- stan(model_code=mod, data = ", label, ".stand, iter = ", iter, ", chains=1, chain_id=", i, ", refresh=1,", sep=""),
                 paste("        control = list(max_treedepth = ", max.tree.depth, "), ", sep=""),
                 paste("        sample_file = '", label, "-fit", i, "_sample.csv', diagnostic_file = '", label, "-fit", i, "_diagnostic.csv')", sep=""),
                 paste("save(", label, ".fit", i, ", file=paste(Sys.getenv('WORK'), '/", clust.path, label, "-fit", i, ".RData', sep=''))", sep="")), fileConn)
    close(fileConn)

    fileConn <- file(paste(file.path, label, "-fit", i, ".pbs", sep=""))
    writeLines(c("#!/bin/sh", 
                 "#PBS -l walltime=72:00:00",
                 "#PBS -l select=01:ncpus=1:mem=4500mb",
                 "#PBS -j oe",
                 "",
                 "module load R/3.2.0",
                 "module load intel-suite/2015.3",
                 "",
                 paste("R CMD BATCH --no-restore --no-save $HOME/", clust.path, label, "-fit", i, ".R $WORK/", clust.path, label, "-fit", i, ".Rout", sep=""),
                 "",
                 "qstat -f $PBS_JOBID"), fileConn)
    close(fileConn)

    fileConn <- file(paste(file.path, label, "-fit", i, ".bat", sep=""))
    writeLines(c(## "net use Q: \\\\fi--san02.dide.ic.ac.uk\\homes\\jwe08",
                 ## "call Q:\\cluster-config.bat",
                 "net use /y P: \\\\fi--didenas1-app\\jeff",
                 "call P:\\cluster-config.bat",
                 paste("Rcmd BATCH --no-restore --no-save %HOME%", sub("/", "\\\\", clust.path), label, "-fit", i, ".R %WORK%", sub("/", "\\\\", clust.path), label, "-fit", i, ".Rout", sep="")),
               fileConn)
    close(fileConn)

    fileConn <- file(paste(file.path, "submit-jobs.sh", sep=""), "a")
    writeLines(paste("qsub ", label, "-fit", i, ".pbs", sep=""), fileConn)
    close(fileConn)

    fileConn <- file(paste(file.path, "submit-jobs.bat", sep=""), "a")
    ## writeLines(paste("job submit /scheduler:fi--dideclusthn.dide.ic.ac.uk /jobtemplate:GeneralNodes /jobname:", label, "-fit", i, " /numcores:1 \\\\fi--san02.dide.ic.ac.uk\\homes\\jwe08\\", sub("/", "\\\\", clust.path), label, "-fit", i, ".bat", sep=""), fileConn)
    writeLines(paste("job submit /scheduler:fi--didemrchnb.dide.ic.ac.uk /jobtemplate:GeneralNodes /jobname:", label, "-fit", i, " /numcores:1 \\\\fi--didenas1.dide.ic.ac.uk\\Jeff\\", sub("/", "\\\\", clust.path), label, "-fit", i, ".bat", sep=""), fileConn)
    close(fileConn)
  }

  return(invisible())
}
