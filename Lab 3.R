library(data.table)
library(ggplot2)

#Question 1
theta <- seq(0, 1, length.out = 1000)
post_A <- dbeta(theta, 0.5 + 11, 0.5 + 3)
post_B <- dbeta(theta, 0.5 + 5, 0.5 + 1)
data_A <- data.table(theta = theta, post = post_A, solution = "A")
data_B <- data.table(theta = theta, post = post_B, solution = "B")
data <- rbindlist(list(data_A, data_B))
ggplot(data, aes(x = theta, y = post, group = solution)) + geom_line(aes(color = solution)) + labs(title = "Posterior Distributions", x = "theta", y = "Density")

#Question 2
prob_A <- 1 - pbeta(0.8, 0.5 + 11, 0.5 + 3)
prob_B <- 1 - pbeta(0.8, 0.5 + 5, 0.5 + 1)

#Question 3, P(theta_B - theta_A | y, x)
theta_A <- rbeta(10000, 0.5 + 11, 0.5 + 3)
theta_B <- rbeta(10000, 0.5 + 5, 0.5 + 1)
diff <- theta_B - theta_A
mean(diff > 0)
