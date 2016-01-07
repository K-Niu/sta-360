library(data.table)
library(ggplot2)

#Question 8
exp_loss <- function(c) {
  f_1 <- function(theta) {
    abs(theta - c)*(theta^0.05)*((1-theta)^29)/beta(1.05, 30)
  }
  f_2 <- function(theta) {
    10*abs(theta - c)*(theta^0.05)*((1-theta)^29)/beta(1.05, 30)
  }
  N = 1000
  integrate(f_1, 0, c, subdivisions = N)$value + integrate(f_2, c, 1, subdivisions = N)$value
}

c = seq(0, 0.5, length.out = 1000)
values = as.numeric(lapply(c, exp_loss))
data <- data.table(c = c, values = values)
ggplot(data, aes(x = c, y = values)) + geom_line() + labs(title = "Posterior expected loss for disease prevalence example", x = "c", y = "Posterior expected loss")

#In-class problem
##Beta(1, 299) prior for theta, Beta(1, 1) prior for beta
beta <- rbeta(1000, 1 + 2, 1 + 28)
theta <- rbeta(1000, 1 + 0, 299 + 30)
diff <- beta - theta
data1 <- data.table(diff = diff, prior = "θ ~ Beta(1, 299), β ~ Beta(1, 1)")
ggplot(data1, aes(x = diff)) + geom_density() + labs(title = "Density of β - θ | y, x", x = "β - θ", y = "Density") #Density of β - θ | y, x
sum(diff > 0)/1000 #Estimate of Pr(β > θ | y, x)
diff = sort(diff)
c(diff[25], diff[975]) #95% credible interval for β - θ

#Beta(1, 1) prior for theta, Beta(1, 1) prior for beta
beta <- rbeta(1000, 1 + 2, 1 + 28)
theta <- rbeta(1000, 1 + 0, 1 + 30)
diff <- beta - theta
data2 <- data.table(diff = diff, prior = "θ ~ Beta(1, 1), β ~ Beta(1, 1)")
ggplot(data2, aes(x = diff)) + geom_density() + labs(title = "Density of β - θ | y, x", x = "β - θ", y = "Density")
sum(diff > 0)/1000
diff = sort(diff)
c(diff[25], diff[975])

#Beta(1, 2499) prior for theta, Beta(1, 1) prior for beta
beta <- rbeta(1000, 1 + 2, 1 + 28)
theta <- rbeta(1000, 1 + 0, 2499 + 30)
diff <- beta - theta
data3 <- data.table(diff = diff, prior = "θ ~ Beta(1, 2499), β ~ Beta(1, 1)")
ggplot(data3, aes(x = diff)) + geom_density() + labs(title = "Density of β - θ | y, x", x = "β - θ", y = "Density")
sum(diff > 0)/1000
diff = sort(diff)
c(diff[25], diff[975])

#Beta(2, 1) prior for theta, Beta(1, 1) prior for beta
beta <- rbeta(1000, 1 + 2, 1 + 28)
theta <- rbeta(1000, 2 + 0, 1 + 30)
diff <- beta - theta
data4 <- data.table(diff = diff, prior = "θ ~ Beta(2, 1), β ~ Beta(1, 1)")
ggplot(data4, aes(x = diff)) + geom_density() + labs(title = "Density of β - θ | y, x", x = "β - θ", y = "Density")
sum(diff > 0)/1000
diff = sort(diff)
c(diff[25], diff[975])

#Compare posteriors
data <- rbindlist(list(data1, data2, data3, data4))
ggplot(data, aes(x = diff, group = prior)) + geom_density(aes(group = prior, color = prior, fill = prior), alpha = 0.3) + labs(title = "Posterior densities for different priors", x = "β - θ", y = "Density")
