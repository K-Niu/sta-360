library(data.table)
library(ggplot2)

#Samples the marginal distribution of X
sample_marginal <- function(nu, num_samples) {
  tao_sq = rgamma(num_samples, nu/2, rate = nu/2)
  normal_with_tao <- function(tao_sq) {
    rnorm(1, 0, sqrt(1/tao_sq))
  }
  x = as.numeric(lapply(tao_sq, normal_with_tao))
}

#Plotting a sample of 10000 from the marginal distribution of X
data <- data.table(x = sample_marginal(1, 10000))
ggplot(data, aes(x = x)) + geom_histogram() + scale_x_continuous(limits = c(-4, 4)) + labs(title = "10000 samples from the marginal distribution of X")

#Test whether observed distribution is equal to a t distribution with 1 df
ks.test(data$x, "pt", 1)

#Sampling and performing KS test 1000 times, using 100 draws from the marginal distribution of X
pvalues = data.table(p = numeric())
for(i in 1:1000) {
  x = sample_marginal(1, 100)
  ks = ks.test(x, "pt", 1)$p.value
  pvalues = rbindlist(list(pvalues, data.table(p = ks)))
}
ggplot(pvalues, aes(x = p)) + geom_histogram() + labs(title = "Distribution of p-values from Kolmogorov-Smirnov test")
