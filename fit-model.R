setwd("~/Documents/Research/ALPHA/life-expectancy/alpha-le-stan/")

source("le-model-functions.R")
load("../../alpha-incidence/episode-data-formatted_20150204.RData")
sm.incmat.pen <- stan_model("le-mod_incmat-pen.stan")



umkhanm.dat <- prepare.data(subset(datw, site=="uMkhanyakude" & sex=="Male"),
                            natmxstart.time=2000.5, artstart.time=2005.5)
umkhanf.dat <- prepare.data(subset(datw, site=="uMkhanyakude" & sex=="Female"),
                            natmxstart.time=2000.5, artstart.time=2005.5)

masakam.dat <- prepare.data(subset(datw, site=="Masaka" & sex=="Male"),
                            natmxstart.time=1990.0, artstart.time=2005.5)
masakaf.dat <- prepare.data(subset(datw, site=="Masaka" & sex=="Female"),
                            natmxstart.time=1990.0, artstart.time=2005.5)

umkhanm.mle <- optimizing(sm.incmat.pen, data=umkhanm.dat, hessian=TRUE, iter=2000)
oumkhanf.mle <- optimizing(sm.incmat.pen, data=umkhanf.dat, hessian=TRUE, iter=2000)

masakam.mle <- optimizing(sm.incmat.pen, data=masakam.dat, hessian=TRUE, iter=2000)



masakaf.dat <- prepare.data(subset(datw, site=="Masaka" & sex=="Female"),
                            natmxstart.time=1990.0, artstart.time=2005.5)

masakaf.mle <- optimizing(sm.incmat.pen, data=masakaf.dat, iter=2000, hessian=TRUE)
