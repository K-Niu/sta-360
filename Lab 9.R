library(data.table)
library(truncnorm)

y <- read.csv('Lab 9 - data.txt', sep = ' ')
y = data.table(y)
y[, weekend := ifelse(day %in% c('Saturday', 'Sunday'), 1, 0)]
y.weekend <- y[weekend == 1,]$pm
y.weekday <- y[weekend == 0,]$pm
n.weekend <- length(y.weekend)
n.weekday <- length(y.weekday)

#Priors
m <- 1
v <- 100
mu.ep <- 1
sigma.sq.ep <- 100
a <- 0.001
b <- 0.001
#Gibbs sampling
mu2.Gibbs <- 1
sigma1.sq.inv.Gibbs <- 1
sigma2.sq.inv.Gibbs <- 1
epsilon.Gibbs <- 1
for(i in 2:10010) {
  mu2.Gibbs[i] = rnorm(1, (1/(sigma1.sq.inv.Gibbs[i-1]*n.weekday + sigma2.sq.inv.Gibbs[i-1]*n.weekend + 1/v))*(sigma1.sq.inv.Gibbs[i-1]*sum(log(y.weekday)) - sigma1.sq.inv.Gibbs[i-1]*epsilon.Gibbs[i-1]*n.weekday + sigma2.sq.inv.Gibbs[i-1]*sum(log(y.weekend)) + m/v),
                       sqrt(1/(sigma1.sq.inv.Gibbs[i-1]*n.weekday + sigma2.sq.inv.Gibbs[i-1]*n.weekend + 1/v)))
  sigma1.sq.inv.Gibbs[i] = rgamma(1, a + 0.5*n.weekday, rate = b + 0.5*sum((log(y.weekday) - (mu2.Gibbs[i] + epsilon.Gibbs[i-1]))^2))
  sigma2.sq.inv.Gibbs[i] = rgamma(1, a + 0.5*n.weekend, rate = b + 0.5*sum((log(y.weekend) - mu2.Gibbs[i])^2))
  epsilon.Gibbs[i] = rtruncnorm(1, a = 0, b = Inf, mean = (1/(sigma1.sq.inv.Gibbs[i]*n.weekday + 1/sigma.sq.ep))*(sigma1.sq.inv.Gibbs[i]*sum(log(y.weekday)) - sigma1.sq.inv.Gibbs[i]*mu2.Gibbs[i]*n.weekday + mu.ep/sigma.sq.ep), sd = sqrt(1/(sigma1.sq.inv.Gibbs[i]*n.weekday + 1/sigma.sq.ep)))
}
#Trace plots
plot(seq(1, 10010), mu2.Gibbs + epsilon.Gibbs, type = 'l', main = 'Trace plot for mu1', xlab = 'iteration', ylab = 'mu1')
plot(seq(1, 10010), mu2.Gibbs, type = 'l', main = 'Trace plot for mu2', xlab = 'iteration', ylab = 'mu2')
plot(seq(1, 10010), 1/sigma1.sq.inv.Gibbs, type = 'l', main = 'Trace plot for sigma1^2', xlab = 'iteration', ylab = 'sigma1^2')
plot(seq(1, 10010), 1/sigma2.sq.inv.Gibbs, type = 'l', main = 'Trace plot for sigma2^2', xlab = 'iteration', ylab = 'sigma2^2')
#Take out burn-ins
mu2.Gibbs = mu2.Gibbs[11:10010]
sigma1.sq.inv.Gibbs = sigma1.sq.inv.Gibbs[11:10010]
sigma2.sq.inv.Gibbs = sigma2.sq.inv.Gibbs[11:10010]
epsilon.Gibbs = epsilon.Gibbs[11:10010]

#Posterior point esimates
mean(mu2.Gibbs + epsilon.Gibbs)
quantile(mu2.Gibbs + epsilon.Gibbs, c(0.025, 0.975))
mean(1/sigma1.sq.inv.Gibbs)
quantile(1/sigma1.sq.inv.Gibbs, c(0.025, 0.975))
mean(mu2.Gibbs)
quantile(mu2.Gibbs, c(0.025, 0.975))
mean(1/sigma2.sq.inv.Gibbs)
quantile(1/sigma2.sq.inv.Gibbs, c(0.025, 0.975))

#Probability mu1 > mu2
mean(epsilon.Gibbs > 0)
#Probability sigma1 > sigma2
mean(sqrt(1/sigma1.sq.inv.Gibbs) > sqrt(1/sigma2.sq.inv.Gibbs))

#Probability pollution on Tue > Saturday
tue <- NULL
sat <- NULL
for(i in 1:10000) {
  tue[i] = rlnorm(1, mu2.Gibbs[i] + epsilon.Gibbs[i], sqrt(1/sigma1.sq.inv.Gibbs[i]))
  sat[i] = rlnorm(1, mu2.Gibbs[i], sqrt(1/sigma2.sq.inv.Gibbs[i]))
}
mean(tue > sat)
