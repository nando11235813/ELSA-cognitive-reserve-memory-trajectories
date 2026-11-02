# Bayesian nonlinear mixed-effects model with a smooth transition to exponential decay
# Self-contained example: SIMULATED data only (no ELSA data).
# Requirements: R, JAGS, and R packages rjags, snow, coda.
# The model uses a constant pre-onset level, transition width dt, and exponential decay.

library(rjags)
library(snow)
library(coda)

# -------------------------------------------------------------------------
# 1. SIMULATED LONGITUDINAL DATA
# -------------------------------------------------------------------------
# Structure of `data` (long format; one row per individual and assessment):
#   id   : integer subject ID, 1, ..., n; repeated across assessments
#   age  : age at assessment (50 to 90 in this example)
#   y    : observed continuous outcome
#   x1   : numerical, time-invariant subject-level covariate
#   x2   : binary (0/1), time-invariant subject-level covariate
# The covariates are repeated across visits of the same subject.
# A real dataset may replace this entire simulation block, provided its
# columns have the same names and id is coded 1,...,n.

set.seed(2026)
n <- 50

# Subject-level covariates and varying number of assessments (4 to 8)
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rbinom(n, size = 1, prob = 0.5)
n_visit <- sample(4:8, size = n, replace = TRUE)

# Subject-level parameters used to SIMULATE the outcome (not estimated here)
b0_true <- rnorm(n, mean = 11 + 0.8*x1 + 1.0*x2, sd = 0.9)
CP_true <- rnorm(n, mean = 68 + 1.0*x1 + 2.0*x2, sd = 2)
b2_true <- rnorm(n, mean = -0.045 - 0.005*x1 - 0.01*x2, sd = 0.008)
dt_true <- 3
sigma_eps_true <- 1.2

data <- do.call(rbind, lapply(seq_len(n), function(j) {
  age_j <- sort(sample(seq(50, 90, by = 1), size = n_visit[j], replace = FALSE))

  # Integral of smoothstep relative-decay rate: 0 before CP, b2 after CP+dt.
  u <- pmin(1, pmax(0, (age_j - CP_true[j])/dt_true))
  cumulative_rate <- b2_true[j] * (
    dt_true*(u^3 - 0.5*u^4) + pmax(age_j - CP_true[j] - dt_true, 0)
  )
  mu_j <- b0_true[j] * exp(cumulative_rate)

  data.frame(id = j,
             age = age_j,
             y = rnorm(length(age_j), mean = mu_j, sd = sigma_eps_true),
             x1 = x1[j],
             x2 = x2[j])
}))
rownames(data) <- NULL

# Illustration: nrow(data) observations, with different visit counts by id.
print(head(data, 10))
print(table(data$id))

# -------------------------------------------------------------------------
# 2. JAGS MODEL: independently distributed regression coefficients
# -------------------------------------------------------------------------
# b0_i: pre-onset level; CP_i: onset; b2_i: post-transition relative rate.
# The same two covariates enter each parameter with separate coefficients.
# Priors preserve the parameterization of the original example.

model_file <- "dem.bug"
writeLines('
model {
  # Half-Cauchy priors for standard deviations (precision = scale^-2)
  s_b0  ~ dt(0, pow(prior_scale[1], -2), 1) T(0, )
  s_b2  ~ dt(0, pow(prior_scale[2], -2), 1) T(0, )
  s_CP  ~ dt(0, pow(prior_scale[3], -2), 1) T(0, )
  s_eps ~ dt(0, pow(prior_scale[4], -2), 1) T(0, )

  # Population-level parameters
  b0 ~ dunif(0, 20)           # level
  b2 ~ dnorm(0, 0.01)         # relative decay rate
  CP ~ dunif(60, 80)          # onset age
  dt ~ dunif(0, 5)            # transition duration

  # Fixed effects: independent priors for each covariate and parameter
  b0_x1 ~ dnorm(0, 0.01)
  b0_x2 ~ dnorm(0, 0.01)
  CP_x1 ~ dnorm(0, 0.01)
  CP_x2 ~ dnorm(0, 0.01)
  b2_x1 ~ dnorm(0, 0.01)
  b2_x2 ~ dnorm(0, 0.01)

  # Individual-level parameters
  for (j in 1:n) {
    b0_i[j] ~ dnorm(b0 + b0_x1*x1[j] + b0_x2*x2[j], pow(s_b0, -2))
    b1_i[j] <- 0
    b2_i[j] ~ dnorm(b2 + b2_x1*x1[j] + b2_x2*x2[j], pow(s_b2, -2))
    CP_i[j] ~ dnorm(CP + CP_x1*x1[j] + CP_x2*x2[j], pow(s_CP, -2))

    # Cubic relative-rate function for the smooth transition.
    # Keep the parameterization used in the original JAGS script.
    d_i[j] <- 2*(b1_i[j]/b0_i[j] - b2_i[j]) / pow(dt, 3)
    c_i[j] <- -3*(2*CP_i[j]+dt)*(b1_i[j]/b0_i[j]-b2_i[j]) / pow(dt, 3)
    b_i[j] <- 6*CP_i[j]*(CP_i[j]+dt)*(b1_i[j]/b0_i[j]-b2_i[j]) / pow(dt, 3)
    a_i[j] <- b1_i[j]/b0_i[j] -
              (2*CP_i[j]+3*dt)*pow(CP_i[j], 2)*(b1_i[j]/b0_i[j]-b2_i[j]) / pow(dt, 3)
  }

  # Longitudinal observations
  for (i in 1:N) {
    mu[i] <- b0_i[gr[i]] * ifelse(age[i] <= CP_i[gr[i]], 1,
      ifelse(age[i] <= CP_i[gr[i]] + dt,
        exp(d_i[gr[i]]*(pow(age[i],4)-pow(CP_i[gr[i]],4))/4 +
            c_i[gr[i]]*(pow(age[i],3)-pow(CP_i[gr[i]],3))/3 +
            b_i[gr[i]]*(pow(age[i],2)-pow(CP_i[gr[i]],2))/2 +
            a_i[gr[i]]*(age[i]-CP_i[gr[i]])),
        exp(b2_i[gr[i]]*(age[i]-(CP_i[gr[i]]+dt)) +
            d_i[gr[i]]*(pow(CP_i[gr[i]]+dt,4)-pow(CP_i[gr[i]],4))/4 +
            c_i[gr[i]]*(pow(CP_i[gr[i]]+dt,3)-pow(CP_i[gr[i]],3))/3 +
            b_i[gr[i]]*(pow(CP_i[gr[i]]+dt,2)-pow(CP_i[gr[i]],2))/2 +
            a_i[gr[i]]*dt)))
    Y[i] ~ dnorm(mu[i], pow(s_eps, -2))
  }
}
', con = model_file)

# -------------------------------------------------------------------------
# 3. ESTIMATION: four independent chains using snow
# -------------------------------------------------------------------------
# As in the original code, covariates are passed once per individual.
# The data supplied to JAGS must include `prior_scale` (required by priors).

n_chains <- 4
n_iter <- 1000       # saved draws per chain
n_iterBI <- 10000    # burn-in iterations per chain
prior_scale <- c(1, 2, 0.1, 1)  # b0, b2, CP, residual: ORIGINAL example

stopifnot(identical(sort(unique(data$id)), seq_len(n)))
stopifnot(all(vapply(c("x1", "x2"), function(v) {
  all(vapply(split(data[[v]], data$id), function(z) length(unique(z)) == 1L,
             logical(1)))
}, logical(1))))

jags_data <- list(
  N = nrow(data),
  n = n,
  Y = data$y,
  gr = data$id,
  age = data$age,
  x1 = data$x1[!duplicated(data$id)],
  x2 = data$x2[!duplicated(data$id)],
  prior_scale = prior_scale
)

# One worker = one MCMC chain. Worker gets model file and data in its environment.
DEM_snow <- function(j) {
  mod <- jags.model(
    model_file,
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1000 + j),
    data = jags_data,
    n.chains = 1,
    n.adapt = 1000,
    quiet = TRUE
  )
  update(mod, n.iter = n_iterBI, progress.bar = "none")
  coda.samples(
    mod,
    variable.names = c(
      "b0", "CP", "b2", "dt",
      "b0_x1", "b0_x2", "CP_x1", "CP_x2", "b2_x1", "b2_x2",
      "s_b0", "s_CP", "s_b2", "s_eps"
    ),
    n.iter = n_iter,
    progress.bar = "none"
  )
}

clus <- makeCluster(n_chains, type = "SOCK")
result <- tryCatch({
  clusterEvalQ(clus, library(rjags))
  clusterExport(clus, c("jags_data", "n_iter", "n_iterBI", "model_file"))
  clusterApply(clus, seq_len(n_chains), DEM_snow)
}, finally = {
  stopCluster(clus)
})

# Each worker returns a one-chain mcmc.list
samples <- do.call(coda::mcmc.list, lapply(result, function(z) z[[1]]))

print(summary(samples))
print(coda::gelman.diag(samples, multivariate = FALSE))
print(coda::effectiveSize(samples))

# Optional: save estimates, but NEVER upload restricted source datasets.
# saveRDS(samples, "posterior_samples_example.rds")
