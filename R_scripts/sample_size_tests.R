
#Sample size sims



# --- paternity
alpha <- 0.05 # significance level
target_power <- 0.8 # desired power (1 - beta)
n_sims <- 200 # monte carlo replicates, increase for precision
sires <- 17 # number of candidate sires
p0 <- rep(1/sires, sires) # null equal probabilities
# EXAMPLES of alternative true sire probs (replace with your scenario)
#alt_small <- c(0.10, rep((1-0.10)/(sires-1), sires-1)) # one sire slightly elevated
skew = 0.35
alt_moderate <- c(skew, rep((1-skew)/(sires-1), sires-1)) # one sire moderately elevated
alt_dirichlet <- as.numeric(rlang::rep_along(seq_len(sires),1)) # placeholder; replace with explicit probs
n_grid <- seq(10, 50, by=1) # grid of larvae/sample sizes to test
# --- end editable block --- 
power_for_alt <- function(alt,p0,n_grid,n_sims,alpha){set.seed(42) # reproducible
  start_time <- Sys.time() # tic: record start time
  power_vec <- numeric(length(n_grid)) # store power for each sample size
  progress_step <- max(1, ceiling(length(n_grid)/20)) # print ~20 updates over the run
  for(i in seq_along(n_grid)){m <- n_grid[i] # larvae per brood to test
  rej <- 0L # count rejections
  for(sim in seq_len(n_sims)){x <- rmultinom(1,m,prob=alt) # draw counts
  test <- suppressWarnings(chisq.test(x,p=p0,simulate.p.value=FALSE)) # Pearson chi-square
  if(is.na(test$p.value)) next # skip if something odd
  if(test$p.value < alpha) rej <- rej + 1L
  }
  power_vec[i] <- rej/n_sims # store power for this m
  if(i %% progress_step == 0){elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs")) # elapsed secs
  message(sprintf("progress: %d/%d (m=%d) — power=%.3f — elapsed=%.0f s", i, length(n_grid), m, power_vec[i], elapsed))} # progress msg
  }
  end_time <- Sys.time() # toc: record end time
  total_secs <- as.numeric(difftime(end_time, start_time, units="secs")) # total seconds
  message(sprintf("done: tested %d sample sizes in %.0f s (≈ %.1f s per grid point)", length(n_grid), total_secs, total_secs/length(n_grid))) # final summary
  data.frame(m = n_grid, power = power_vec) # return dataframe of power vs m
}

find_min_n <- function(power_df,target_power){idx <- which(power_df$power>=target_power)
if(length(idx)==0) return(NA) else return(min(power_df$m[idx]))}
# run for example alternatives (replace alt_* with the actual alt you want)
#res_small <- power_for_alt(alt_small,p0,n_grid,n_sims,alpha) # small skew
res_moderate <- power_for_alt(alt_moderate,p0,n_grid,n_sims,alpha) # moderate skew
# print minima
#cat("min larvae for small skew:",find_min_n(res_small,target_power),"\n") # report
cat("min larvae for moderate skew:",find_min_n(res_moderate,target_power),"\n") # report



# neg poisson -------------------------------------------------------------


# library(glmmTMB)
# set.seed(42)
# simulate_power_nb_glmm <- function(n_pairs, beta_dist, beta_angle, theta, n_sims=200, alpha=0.05) {
#   pvals <- numeric(n_sims)
#   start_time <- Sys.time() # tic: start timer
#   progress_step <- max(1, floor(n_sims / 10)) # print progress every 10% of sims
#   for(i in seq_len(n_sims)) {
#     dist_m_c <- runif(n_pairs, min=-2, max=2) # example predictor range; adjust as needed
#     cos_ang_c <- runif(n_pairs, min=-1, max=1) # example predictor range
#     eta <- beta_dist * dist_m_c + beta_angle * cos_ang_c
#     mu <- exp(eta)
#     y <- rnbinom(n_pairs, mu=mu, size=theta) # generate counts
#     dat <- data.frame(final_count = y, dist_m_c = dist_m_c, cos_ang_c = cos_ang_c)
#     fit <- try(glmmTMB(final_count ~ dist_m_c + cos_ang_c, family=nbinom2(), data=dat), silent=TRUE)
#     if(inherits(fit, "try-error")) next
#     sum_fit <- summary(fit)
#     pval_dist <- sum_fit$coefficients$cond["dist_m_c", "Pr(>|z|)"]
#     pvals[i] <- ifelse(is.na(pval_dist), 1, pval_dist) # store p-value for dist_m_c
#     if(i %% progress_step == 0) {
#       elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
#       message(sprintf("Sim %d/%d done (%.1f%%) — elapsed time: %.0f seconds", i, n_sims, 100 * i/n_sims, elapsed))
#     }
#   }
#   end_time <- Sys.time() # toc: end timer
#   total_time <- as.numeric(difftime(end_time, start_time, units="secs"))
#   message(sprintf("Simulation complete: %d sims in %.0f seconds (%.2f s per sim)", n_sims, total_time, total_time/n_sims))
#   mean(pvals < alpha) # proportion significant = empirical power
# }
# # Example usage:
# theta <- 1.5 # dispersion estimated from your model
# beta_dist <- -0.3 # hypothetical effect size for dist_m_c (negative means counts decrease with distance)
# beta_angle <- 0.2 # hypothetical effect for cos_ang_c
# (power_408 <- simulate_power_nb_glmm(408, beta_dist, beta_angle, theta))
# 
# 
# 
# library(stats)
# set.seed(42)
# 
# simulate_power_binom_glm <- function(n_females, beta_intercept, beta_dist, beta_deg, tot_eggs=100, n_sims=200, alpha=0.05) {
#   pvals <- numeric(n_sims)
#   for (i in seq_len(n_sims)) {
#     dist <- runif(n_females, min=0, max=20) # example distance range in metres
#     deg <- runif(n_females, min=-180, max=180) # example degree range
#     eta <- beta_intercept + beta_dist * dist + beta_deg * deg
#     p <- plogis(eta) # convert linear predictor to probability scale
#     
#     suc <- rbinom(n_females, size=tot_eggs, prob=p) # successes per female
#     dat <- data.frame(suc=suc, fail=tot_eggs - suc, dist=dist, deg=deg)
#     
#     fit <- try(glm(cbind(suc, fail) ~ dist + deg, family=binomial, data=dat), silent=TRUE)
#     if (inherits(fit, "try-error")) next
#     sum_fit <- summary(fit)
#     pval_dist <- coef(sum_fit)[2, "Pr(>|z|)"] # p-value for dist coefficient
#     pvals[i] <- ifelse(is.na(pval_dist), 1, pval_dist)
#   }
#   mean(pvals < alpha) # empirical power to detect effect of dist
# }




# sample size fert success ------------------------------------------------

# You can loop over n_females to find sample size for desired power

simulate_power_binom_glm <- function(n_females, beta_intercept, beta_dist, tot_eggs=200, n_sims=200, alpha=0.05) {
  pvals <- numeric(n_sims)
  fixed_dist <- rep(c(0.7, 3, 10, 30), each = 6) # fixed distances for colonies (6 each)
  for (i in seq_len(n_sims)) {
    dist <- fixed_dist[seq_len(n_females)] # subset in case n_females < 24
    eta <- beta_intercept + beta_dist * dist
    p <- plogis(eta)
    suc <- rbinom(n_females, size=tot_eggs, prob=p)
    dat <- data.frame(suc=suc, fail=tot_eggs - suc, dist=dist)
    fit <- try(glm(cbind(suc, fail) ~ dist, family=binomial, data=dat), silent=TRUE)
    if (inherits(fit, "try-error")) next
    sum_fit <- summary(fit)
    pval_dist <- coef(sum_fit)[2, "Pr(>|z|)"]
    pvals[i] <- ifelse(is.na(pval_dist), 1, pval_dist)
  }
  mean(pvals < alpha)
}


# Parameters for 80% at 0 m and 72% at 30 m:


# Calculate intercept and slope for 80% at 0 m and 70% at 30 m

p0 <- 0.80
p30 <- 0.70
beta_intercept <- qlogis(p0) # logit of initial fertilisation probability
beta_dist <- (qlogis(p30) - qlogis(p0)) / 30 # slope for 10% drop over 30m

# Print to check values
cat("beta_intercept =", round(beta_intercept, 3), "\n") # ~1.386
cat("beta_dist =", round(beta_dist, 4), "\n") # ~-0.018

power <- simulate_power_binom_glm(19, beta_intercept, beta_dist)
cat("Power with X females:", power, "\n")


