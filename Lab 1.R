n <- 30
p <- 0.7
i <- 10000

data <- rbinom(i, n, p)
freq_interval_lengths <- NULL
unif_interval_lengths <- NULL
beta_interval_lengths <- NULL
freq_interval_captures <- NULL
unif_interval_captures <- NULL
beta_interval_captures <- NULL

for(i in 1:length(data)) {
  successes <- data[i]
  x_bar <- successes/n
  interval_freq <- c(x_bar - 1.96*sqrt((x_bar*(1 - x_bar))/n), x_bar + 1.96*sqrt((x_bar*(1 - x_bar))/n))
  interval_prior_unif <- qbeta(c(0.025, 0.975), 1 + successes, 1 + n - successes)
  interval_prior_beta <- qbeta(c(0.025, 0.975), 5 + successes, 2 + n - successes)
  freq_interval_lengths = cbind(freq_interval_lengths, interval_freq[2] - interval_freq[1])
  unif_interval_lengths = cbind(unif_interval_lengths, interval_prior_unif[2] - interval_prior_unif[1])
  beta_interval_lengths = cbind(beta_interval_lengths, interval_prior_beta[2] - interval_prior_beta[1])
  freq_interval_captures = cbind(freq_interval_captures, p > interval_freq[1] & p < interval_freq[2])
  unif_interval_captures = cbind(unif_interval_captures, p > interval_prior_unif[1] & p < interval_prior_unif[2])
  beta_interval_captures = cbind(beta_interval_captures, p > interval_prior_beta[1] & p < interval_prior_beta[2])
}

freq_coverage <- mean(freq_interval_captures)
unif_coverage <- mean(unif_interval_captures)
beta_coverage <- mean(beta_interval_captures)
average_freq_length <- mean(freq_interval_lengths)
average_unif_length <- mean(unif_interval_lengths)
average_beta_length <- mean(beta_interval_lengths)

