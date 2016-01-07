x <- read.csv('Lab 8 - data.txt')$pm

#Priors
mu.0 <- 0
sigma.sq.0 <- 100
a <- 0.001
b <- 0.001

#Sampling
mu.Gibbs <- 0
sigma.sq.inv.Gibbs <- 1
n <- length(x)
for(i in 2:10010) {
  mu.Gibbs[i] = rnorm(1, (1/(n*sigma.sq.inv.Gibbs[i-1] + 1/sigma.sq.0))*(sigma.sq.inv.Gibbs[i-1]*sum(log(x)) + (1/sigma.sq.0)*mu.0), sqrt(1/(n*sigma.sq.inv.Gibbs[i-1] + 1/sigma.sq.0)))
  sigma.sq.inv.Gibbs[i] = rgamma(1, a + n/2, b + 0.5*sum((log(x) - mu.Gibbs[i])^2))
}

#Take out burn-ins
plot(seq(1, 10010), mu.Gibbs, type = 'l', main = 'Trace plot for mu', xlab = 'iteration')
plot(seq(1, 10010), sigma.sq.inv.Gibbs, type = 'l', main = 'Trace plot for sigma^(-2)', xlab = 'iteration')
mu.Gibbs = mu.Gibbs[11:10010]
sigma.sq.inv.Gibbs = sigma.sq.inv.Gibbs[11:10010]

#95% CI for mean and variance of x
x.mean <- sapply(seq(1, 10000), function(i, mu, sigma.sq.inv) {
  return(exp(mu[i] + 1/(2*sigma.sq.inv[i])))
}, mu = mu.Gibbs, sigma.sq.inv = sigma.sq.inv.Gibbs)
x.var <- sapply(seq(1, 10000), function(i, mu, sigma.sq.inv) {
  (exp(1/sigma.sq.inv[i]) - 1)*exp(2*mu[i] + 1/sigma.sq.inv[i])
}, mu = mu.Gibbs, sigma.sq.inv = sigma.sq.inv.Gibbs)
x.mean.interval <- quantile(x.mean, c(0.025, 0.975))
x.var.interval <- quantile(x.var, c(0.025, 0.975))
