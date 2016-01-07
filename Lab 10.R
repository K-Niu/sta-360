library(truncnorm)

y <- read.csv('Lab 10 - data.txt')$pm
n = length(y)

#Priors
m <- 1
v <- 100
mu.ep <- 1
sigma.sq.ep <- 100
a <- 0.001
b <- 0.001
#Gibbs sampling
c.Gibbs <- list(rep(0, n))
p.Gibbs <- 5/7
mu2.Gibbs <- 1
epsilon.Gibbs <- 1
sigma1.sq.inv.Gibbs <- 1
sigma2.sq.inv.Gibbs <- 1
for(k in 2:10010) {
  c.Gibbs[[k]] = sapply(seq(1, n), function(i) {
    prob = (p.Gibbs[k-1]*dlnorm(y[i], mu2.Gibbs[k-1] + epsilon.Gibbs[k-1], sqrt(1/sigma1.sq.inv.Gibbs[k-1])))/(p.Gibbs[k-1]*dlnorm(y[i], mu2.Gibbs[k-1] + epsilon.Gibbs[k-1], sqrt(1/sigma1.sq.inv.Gibbs[k-1])) + (1 - p.Gibbs[k-1])*dlnorm(y[i], mu2.Gibbs[k-1], sqrt(1/sigma2.sq.inv.Gibbs[k-1])))
    return(rbinom(1, 1, prob))
  })
  p.Gibbs[k] = rbeta(1, 5 + sum(c.Gibbs[[k]]), 2 + sum(1 - c.Gibbs[[k]]))
  mu2.Gibbs[k] = rnorm(1, (1/(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]) + sigma2.sq.inv.Gibbs[k-1]*sum(1 - c.Gibbs[[k]]) + 1/v))*(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]*log(y)) - sigma1.sq.inv.Gibbs[k-1]*epsilon.Gibbs[k-1]*sum(c.Gibbs[[k]]) + sigma2.sq.inv.Gibbs[k-1]*sum((1 - c.Gibbs[[k]])*log(y)) + m/v), sqrt(1/(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]) + sigma2.sq.inv.Gibbs[k-1]*sum(1 - c.Gibbs[[k]]) + 1/v)))
  epsilon.Gibbs[k] = rtruncnorm(1, a = 0, b = Inf, mean = (1/(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]) + 1/sigma.sq.ep))*(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]*log(y)) - sigma1.sq.inv.Gibbs[k-1]*mu2.Gibbs[k]*sum(c.Gibbs[[k]]) + mu.ep/sigma.sq.ep), sd = sqrt(1/(sigma1.sq.inv.Gibbs[k-1]*sum(c.Gibbs[[k]]) + 1/sigma.sq.ep)))
  sigma1.sq.inv.Gibbs[k] = rgamma(1, a + 0.5*sum(c.Gibbs[[k]]), b + 0.5*sum(c.Gibbs[[k]]*(log(y) - (mu2.Gibbs[k] + epsilon.Gibbs[k]))^2))
  sigma2.sq.inv.Gibbs[k] = rgamma(1, a + 0.5*sum(1 - c.Gibbs[[k]]), b + 0.5*sum((1 - c.Gibbs[[k]])*(log(y) - mu2.Gibbs[k])^2))
}
#Trace plot
plot(seq(1, 10010), mu2.Gibbs + epsilon.Gibbs, type = 'l', main = "Trace plot for mu1", xlab = 'iteration', ylab = 'mu1')
plot(seq(1, 10010), mu2.Gibbs, type = 'l', main = 'Trace plot for mu2', xlab = 'iteration', ylab = 'mu2')
plot(seq(1, 10010), sigma1.sq.inv.Gibbs, type = 'l', main = 'Trace plot for sigma1^(-2)', xlab = 'iteration', ylab = 'sigma1^(-2)')
plot(seq(1, 10010), sigma2.sq.inv.Gibbs, type = 'l', main = 'Trace plot for sigma2^(-2)', xlab = 'iteration', ylab = 'sigma2^(-2)')

#Take out burn-ins
mu2.Gibbs = mu2.Gibbs[11:10010]
epsilon.Gibbs = epsilon.Gibbs[11:10010]
sigma1.sq.inv.Gibbs = sigma1.sq.inv.Gibbs[11:10010]
sigma2.sq.inv.Gibbs = sigma2.sq.inv.Gibbs[11:10010]

#mu1 vs mu2
plot(mu2.Gibbs + epsilon.Gibbs, mu2.Gibbs, main = 'mu2 vs. mu1', xlab = 'mu1', ylab = 'mu2')

#Number of weekdays
c <- sapply(c.Gibbs, sum)
mean(c)
quantile(c, c(0.025, 0.975))
mean(c/n > 5/7)

#Posterior summaries
mean(mu2.Gibbs + epsilon.Gibbs)
quantile(mu2.Gibbs + epsilon.Gibbs, c(0.025, 0.975))
mean(mu2.Gibbs)
quantile(mu2.Gibbs, c(0.025, 0.975))
1/mean(sigma1.sq.inv.Gibbs)
quantile(1/sigma1.sq.inv.Gibbs, c(0.025, 0.975))
1/mean(sigma2.sq.inv.Gibbs)
quantile(1/sigma2.sq.inv.Gibbs, c(0.025, 0.975))
