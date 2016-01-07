###Homework 1###
library(ggplot2)
library(data.table)

#Question 2
theta <- seq(0, 0.2, 0.001)
prior <- dgamma(theta, 0.1, rate = 1.0)
x <- c(20.9, 69.7, 3.6, 21.8, 21.4, 0.4, 6.7, 10.0)
posterior <- dgamma(theta, 0.1 + length(x), rate = 1.0 + sum(x))
data.prior <- data.table(theta = theta, value = prior, dist = "prior")
data.post <- data.table(theta = theta, value = posterior, dist = "post")
data <- rbindlist(list(data.prior, data.post))

ggplot(data, aes(x = theta, y = value, group = dist)) + geom_line(aes(color = dist)) + labs(title = "Prior and Posterior Distributions")
