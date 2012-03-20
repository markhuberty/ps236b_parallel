###########################
## Demonstration code for 21 March 2012
## PS236b Parallel computation demonstration
## Author: Mark Huberty
## Date: 20 March 2012
###########################

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

      if(n.count == 0 | n.count == 1)
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

    if(i == 0 | i == 1)
      {
        
        return(i)
        
      }else{
        
        out <- fib.seq(N)
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


## Let's try something more pedestrian: bootstrapping

mean.boot.fun <- function(data.vec){

  sample.vec <- sample(1:length(data.vec), replace=TRUE)
  out <- mean(data.vec[sample.vec])
  return(out)

}

boot.seq <- function(N, data.in){

  boots <- sapply(1:N, function(x){

    mean.boot.fun(data.in)

  })

  return(c(mean(boots), quantile(boots, c(0.025, 0.975))))
}

boot.foreach.par <- function(N, data.in){

  boots <- foreach(i=1:N, .combine=c) %dopar% {

    mean.boot.fun(data.in)

  }

  return(c(mean(boots), quantile(boots, c(0.025, 0.975))))
}

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

## A practical example:
## Survey of Consumer Finances 2007
## 22k rows, 5.8k columns

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

## Bootstrap the betas for the regression of
## age on income

lm.fun <- function(x, y){

  sample.vec <- sample(1:length(x), length(x), replace=TRUE)
  x.temp = x[sample.vec]
  y.temp = y[sample.vec]
  lm.temp <- lm(y ~ x)
  lm.beta <- lm$coefficients[2]
  return(lm.beta)
  
}

boot.seq <- function(N, x, y, p.quantiles=c(0.025, 0.975)){
  stopifnot(length(x) == length(y))
  
  out <- sapply(1:N, function(z){

    lm.fun(x, y)
        
  })

  return(mean(out), quantile(out, p.quantiles))

}

boot.foreach.par <-
  function(N, x, y, p.quantiles=c(0.025, 0.975)){

    out <- foreach(i=1:N, .combine=c) %dopar% {

      lm.fun(x, y)

    }

    return(mean(out), quantile(out, p.quantiles))

  }

boot.snow.par <-
  function(N, x, y, p.quantiles=c(0.025, 0.975), cluster){

    out <- parSapply(cluster, 1:N, function(z){

      lm.fun(x, y)

    })

    return(mean(out), quantile(out, p.quantiles))

  }

## Time it:
clusterSetupRNG(cl)
clusterExport(cl, c("lm.fun", "scf"))
time.lmboot.seq <- system.time(boot.seq(1000, scf$age,
                                        scf$inc.normal
                                        )
                               )
time.lmboot.foreach <- system.time(boot.foreach.par(1000,
                                                    scf$age,
                                                    scf$inc.normal
                                                    )
                                   )
time.lmboot.snow <- system.time(boot.snow.par(1000, scf$age,
                                              scf$inc.normal,
                                              cluster=cl
                                              )
                                )
print(time.lmboot.seq)
print(time.lmboot.foreach)
print(time.lmboot.snow)
