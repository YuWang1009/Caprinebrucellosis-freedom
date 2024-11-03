library(tidyverse)
library(R2jags)

# first install jags before running model, available at https://mcmc-jags.sourceforge.io/

rm(list = ls())
gc()

# the below code is adopted from Wood's publication
# Wood, C. et al. (2021). Prevalence and spatial distribution of Coxiella burnetii
# seropositivity in northern Australian beef cattle adjusted for diagnostic test uncertainty.
# Preventive Veterinary Medicine, 189, 105282.

dat <- read.csv("caprinebrucellosis_2021.csv", header = T)

model <- "model{
  for (i in 1:82) {
    y[i] ~ dbin(ap[i], n[i])
    ap[i] <- se * pi[i] + (1 - pi[i]) * (1 - sp)
    logit(pi[i]) <- b + u[county[i]]
  }
  # Priors
  b ~ dnorm(0,0.1)

  for (k in 1:82) {
    u[county[k]] ~ dnorm(w[region[k]], 1 / (sigmau^2))
  }

  for (j in 1:16) {
    w[j] ~ dnorm(muw, 1 / (sigmaw^2))
  }

  seRBT ~ dbeta(14.62, 4.36)T(1-spRBT,)
  spRBT ~ dbeta(5.93, 1.02)
  seSAT ~ dbeta(34.46, 14.01)T(1-spSAT,)
  spSAT ~ dbeta(5.98, 1.04)

  cov_se ~ dunif(0, min(seRBT * (1 - seSAT), seSAT * (1 - seRBT)))
  cov_sp ~ dunif(0, min(spRBT * (1 - spSAT), spSAT * (1 - spRBT)))

  #cov_se ~ dunif(max(-(1 - seRBT) * (1 - seSAT), -seRBT * seSAT), min(seRBT * (1 - seSAT), seSAT * (1 - seRBT)))
  #cov_sp ~ dunif(max(-(1 - spRBT) * (1 - spSAT), -spRBT * spSAT), min(spRBT * (1 - spSAT), spSAT * (1 - spRBT)))

  se <- seRBT * seSAT + cov_se
  sp <- 1 - ((1 - spRBT) * (1 - spSAT) + cov_sp)

  muw <- 0
  sigmau ~ dunif(0, 1)
  sigmaw ~ dunif(0, 1)

  # Calculated outputs
  for (i in 1:16) {
     for (j in 1:82) {
  # prediction for a county in each region
  xb[i,j] <- b + w[i] + u[j]
  # predicted prevalence for a county in each region
  xb_prev_county[i,j] <- exp(xb[i,j]) / (1 + exp(xb[i,j]))
  DisFreePvalue_county[i,j] <- 1 - step(xb_prev_county[i,j] - 0.001)
  }
  xc[i] <- b + w[i]
  # predicted prevalence for a region
  xc_prev_region[i] <- exp(xc[i])/(1+exp(xc[i]))
  DisFreePvalue_region[i] <- 1 - step(xc_prev_region[i] - 0.001)
  }
}"


# Parameters to be monitored
parameters <- c("se", "sp", "xb_prev_county", "DisFreePvalue_county",
                "xc_prev_region","DisFreePvalue_region","seRBT","spRBT","seSAT","spSAT","cov_se","cov_sp")

# Initials
inits1 <- list("seRBT" = 0.8, "spRBT" = 0.9, "seSAT" = 0.8, "spSAT" = 0.9)
inits2 <- list("seRBT" = 0.7, "spRBT" = 0.99, "seSAT" = 0.7, "spSAT" = 0.99)
inits3 <- list("seRBT" = 0.9, "spRBT" = 0.95, "seSAT" = 0.9, "spSAT" = 0.95)
inits <- list(inits1, inits2, inits3)

data_jags <- list(
  y = dat$n_pos,
  n = dat$n_sample,
  region = dat$region,
  county = dat$county
)

# Fitting the model using R JAGS
Sys.setenv(JAGS_HOME="D:/tools/JAGS/JAGS-4.3.1")
library(R2jags)

start <- Sys.time()
bayes.mod.fit <- jags(
  data = data_jags,
  inits = inits,
  parameters.to.save = parameters,
  n.chains = 3,
  n.iter = 250000,
  n.burnin = 50000,
  model.file = textConnection(model)
)

end <- Sys.time()
(end-start)

bayes.mod.fit

traceplot(bayes.mod.fit)
