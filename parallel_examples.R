###########################
## Demonstration code for 21 March 2012
## PS236b Parallel computation demonstration
## Author: Mark Huberty
## Date: 20 March 2012
###########################
setwd("~/Documents/Technology/Documentation/parallel_slides")
## Parallel frameworks
library(doMC)
library(doMPI)
library(doSNOW)
library(snow)
## Parallel looping constructs
library(foreach)
## Useful Hadley Wickham packages
library(reshape)
library(plyr)
library(ggplot2)

## Fibonnaci sequence function

fib.seq <- function(N){

  n.count <- 0
  fib.vec <- c()
  while(n.count < N)
    {

      if(n.count < 2)
        {

          fib.vec <- append(fib.vec, n.count)
          
        }else{

          this.iter <- fib.vec[n.count - 1] + fib.vec[n.count - 2]
          fib.vec <- append(fib.vec, this.iter)
        }

      n.count <- n.count + 1

    }

  return(fib.vec)

}

## Parallel fibonnaci (dumb)
fib.par <- function(N){

  n.vec <- 1:N
  fib.vec <- foreach(i=n.vec, .combine=c) %dopar% {

    if(i < 2)
      {
        
        return(i)
        
      }else{
        
        out <- fib.seq(i)
        return(out[length(out)])
        
      }    
  }

  return(fib.vec)

}


## Time the sequential fib:
time.seq <- system.time(fib.seq(500))

cl <- makeSOCKcluster(c("localhost", "localhost"))
registerDoSNOW(cl)
clusterExport(cl, "fib.seq")
time.par <- system.time(fib.par(500))

print(time.seq)
print(time.par)


###################################################
## Bootstrapping example

## Define the core bootstrap function
mean.boot.fun <- function(data.vec){

  sample.vec <- sample(1:length(data.vec), replace=TRUE)
  out <- mean(data.vec[sample.vec])
  return(out)

}

## Define the sequential bootstrapper
boot.seq <- function(N, data.in){

  boots <- sapply(1:N, function(x){

    mean.boot.fun(data.in)

  })

  return(c(mean(boots), quantile(boots, c(0.025, 0.975))))
}

## Define a parallel bootstrapper based on foreach()
boot.foreach.par <- function(N, data.in){

  boots <- foreach(i=1:N, .combine=c) %dopar% {

    mean.boot.fun(data.in)

  }

  return(c(mean(boots), quantile(boots, c(0.025, 0.975))))
}

## Define a parallel bootstrapper based on snow functions
boot.snow.par <- function(N, data.in, cluster){

  boots <- parSapply(cluster, 1:N, function(x){

    mean.boot.fun(data.in)
    
  })

  return(c(mean(boots), quantile(boots, c(0.025, 0.975))))
                     
}

## Time it:
data.in <- 1:10000
clusterExport(cl, c("mean.boot.fun", "data.in"))
boot.time.seq <- system.time(boot.seq(10000, data.in))
boot.time.foreach <- system.time(boot.foreach.par(10000, data.in))
boot.time.snow <- system.time(boot.snow.par(10000, data.in, cl))

## And print some timings
print(boot.time.seq)
print(boot.time.foreach)
print(boot.time.snow)

## Why so slow? How long does the mean.boot.fun function take?
mbf.time <- sapply(1:1000, function(x){

  system.time(mean.boot.fun(data.in))[2]

})

## The timing is very very short:
mean(mbf.time)

## Short answer:
## foreach() is convenient but imposes overhead
## parLapply() requires more work but is faster
## foreach() is a good choice when the internal function is
## time-intensive (less communication, more doing)

########################################################
## Bootstrapping example 2

## A practical example:
## Survey of Consumer Finances 2007
## 22k rows, 5.8k columns
library(foreign)
scf <- read.dta("scf2007short.dta") ## This is a big file and takes a bit to load



## Take only the columns I need
attach(scf)
scf.subset <- data.frame(x8022, # age
                         x6670, # empl status
                         x7402, # industry of empl
                         x5702, # yearly inc
                         x7650, # income usual?
                         x7362, # more normal income?
                         x5904, # college?
                         x6306, # govt health insurance
                         x6329, # pvt / empl health insurance
                         x7186, # Save for future expenses?
                         x3023  # Savings adequate for retirement?
                         )
detach(scf)
rm(scf)
gc()

## Shorten the object name for convenience and throw out the old one.
scf <- scf.subset
rm(scf.subset)
gc()
names(scf) <- c("age", "empl.status",
                "ind.empl", "inc.yr",
                "inc.usual",
                "inc.normal",
                "college",
                "govt.health.insur",
                "pvt.health.insur",
                "save",
                "save.adq.retire"
                )
scf <- scf[scf$inc.yr > 0,]

## Bootstrap the betas for the regression of
## age on income

lm.fun <- function(x, y){

  sample.vec <- sample(1:length(x), length(x), replace=TRUE)
  x.temp = x[sample.vec]
  y.temp = y[sample.vec]
  lm.temp <- lm(y.temp ~ x.temp)
  lm.beta <- lm.temp$coefficients[2]
  return(lm.beta)
  
}

boot.seq <- function(N, x, y, p.quantiles=c(0.025, 0.975)){
  stopifnot(length(x) == length(y))
  
  out <- sapply(1:N, function(z){

    lm.fun(x, y)
        
  })

  return(c(mean(out), quantile(out, p.quantiles)))

}

boot.foreach.par <-
  function(N, x, y, p.quantiles=c(0.025, 0.975)){

    out <- foreach(i=1:N, .combine=c) %dopar% {

      lm.fun(x, y)

    }

    return(c(mean(out), quantile(out, p.quantiles)))

  }

boot.snow.par <-
  function(N, x, y, p.quantiles=c(0.025, 0.975), cluster){

    out <- parSapply(cluster, 1:N, function(z){

      lm.fun(x, y)

    })

    return(c(mean(out), quantile(out, p.quantiles)))

  }

## Plot the data:
plot(log10(scf$inc.yr) ~ scf$age)

## Time it:
clusterSetupRNG(cl)
clusterExport(cl, c("lm.fun", "scf"))

time.lmboot.seq <- system.time(boot.seq(1000,
                                        scf$age,
                                        log10(scf$inc.yr)
                                        )
                               )

time.lmboot.foreach <- system.time(boot.foreach.par(1000,
                                                    scf$age,
                                                    log10(scf$inc.yr)
                                                    )
                                   )

time.lmboot.snow <- system.time(boot.snow.par(1000,
                                              scf$age,
                                              log10(scf$inc.yr),
                                              cluster=cl
                                              )
                                )

print(time.lmboot.seq)
print(time.lmboot.foreach)
print(time.lmboot.snow)


## Implicit parallelism
## NOTE: pnmath requires gcc-4.2+. Most Macs
## only come w/ the Apple xCode version 4.0.1 or so (maybe
## later for Snow Leopard / Lion?) This code works on my
## linux box...
## Taken from  http://rdav.nics.tennessee.edu/system/files/tgcc-r-handout-2011-07-06.pdf
## v1 <- runif(1000)
## v2 <- runif(100000000)
## system.time(qtukey(v1,2,3))
## system.time(exp(v2))
## system.time(sqrt(v2))

## library(pnmath)
## system.time(qtukey(v1,2,3))
## system.time(exp(v2))
## system.time(sqrt(v2))

## Timings on 4 cores:
## > v1 <- runif(1000)
## > v2 <- runif(100000000)

## > system.time(qtukey(v1,2,3))
##    user  system elapsed 
##  19.640   0.000  19.643 
## > system.time(exp(v2))
##    user  system elapsed 
##   4.200   0.290   4.484 
## > system.time(sqrt(v2))
##    user  system elapsed 
##   2.030   0.220   2.248

## > library(pnmath)
## > system.time(qtukey(v1,2,3))
##    user  system elapsed 
##  19.440   0.000   4.996 
## > system.time(exp(v2))
##    user  system elapsed 
##   4.930   0.380   1.329 
## > system.time(sqrt(v2)
## + )
##    user  system elapsed 
##   2.790   0.310   0.785 

            
