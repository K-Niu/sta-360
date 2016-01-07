library(data.table)
library(ggplot2)

###Question 4
##Setup
mu <- 0
var <- 1
n <- 100000
y <- rnorm(n, mu, sqrt(var))

##Sampling from marginal posterior
#Prior parameters
mu.prior <- 0
k <- 1
a <- 0.001
b <- 0.001
#Joint posterior parameters
mu.post <- (mu.prior + k*sum(y))/(1 + k*n)
a.post <- a + n/2
b.post <- b + (1/2)*sum((y - mean(y))^2) + (n*(mean(y) - mu.prior)^2)/(2*k*(1/k + n))
#Sampling
mu.marginal <- rt(1000, 2*a.post)/sqrt((a.post*(1 + n*k))/(k*b.post)) + mu.post

##Gibbs sampling
#Prior parameters
mu.prior <- 0
var.prior <- 100
tao.prior <- 1/var.prior
a <- 0.001
b <- 0.001
#Sampling
mu.gibbs <- 0
tao.gibbs <- 0
for(i in 2:1010) {
  #mu.p1 = (var.prior/(n*var.prior + 1/tao.gibbs[i-1]))*sum(y) + ((1/tao.gibbs[i-1])/(n*var.prior + 1/tao.gibbs[i-1]))*mu.prior
  #mu.p2 = ((1/tao.gibbs[i-1])*var.prior)/(n*var.prior + 1/tao.gibbs[i-1])
  mu.p1 = (tao.gibbs[i-1]*sum(y) + tao.prior*mu.prior)/(n*tao.gibbs[i-1] + tao.prior)
  mu.p2 = 1/(n*tao.gibbs[i-1] + tao.prior)
  mu.gibbs[i] = rnorm(1, mu.p1, sqrt(mu.p2))
  tao.gibbs[i] = rgamma(1, a + (n/2), b + (1/2)*sum((y - mu.gibbs[i]))^2)
}
#Throw out burn-ins
plot(seq(1, 30), mu.gibbs[1:30], type = "l")
mu.gibbs = mu.gibbs[11:1010]

##Posterior summaries
c(mean(mu.marginal), var(mu.marginal))
c(mean(mu.gibbs), var(mu.gibbs))
mu.marginal = sort(mu.marginal)
mu.gibbs = sort(mu.gibbs)
c(mu.marginal[25], mu.marginal[975])
c(mu.gibbs[25], mu.gibbs[975])

##Visualize
data.marginal <- data.table(mu = mu.marginal, method = "marginal")
data.gibbs <- data.table(mu = mu.gibbs, method = "Gibbs")
data <- rbindlist(list(data.marginal, data.gibbs))
ggplot(data, aes(x = mu, group = method)) + geom_density(aes(fill = method), alpha = 0.4) + scale_x_continuous(limits = c(-0.1, 0.1)) + labs(title = "Comparison of two different choices of priors on posterior density of mu")
ggplot(data.marginal, aes(x = mu)) + geom_density()
ggplot(data.gibbs, aes(x = mu)) + geom_density()
