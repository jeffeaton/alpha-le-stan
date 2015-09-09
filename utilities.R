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

compile.chains <- function(prefix){
  files <- system(paste("ls ", prefix, "*.RData", sep=""), inter=TRUE)
  print(paste(prefix, ": ", length(files), " chains found", sep=""))
  lapply(files, load, environment())
  sflist2stanfit(lapply(grep(prefix, ls(), value=TRUE), get, envir=environment()))
}



## ## Example:
## sourceurl <- "https://raw.github.com/jeffeaton/alpha-le-stan/master/"
## mod <- catfiles("le-functions.stan", "le-data.stan", "le-model.stan", path=sourceurl)
