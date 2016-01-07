library(data.table)
library(ggplot2)
library(MASS)

#Question 1
x.mean <- 0 + c(0.9, 0.1)%*%solve(matrix(c(1, 0.1, 0.1, 1), 2, byrow = TRUE))
x.var <- 1 - c(0.9, 0.1)%*%solve(matrix(c(1, 0.1, 0.1, 1), 2, byrow = TRUE))%*%matrix(c(0.9, 0.1), 2)
y.mean <- 0 + c(0.9, 0.1)%*%solve(matrix(c(1, 0.1, 0.1, 1), 2, byrow = TRUE))
y.var <- 1 - c(0.9, 0.1)%*%solve(matrix(c(1, 0.1, 0.1, 1), 2, byrow = TRUE))%*%matrix(c(0.9, 0.1), 2)
z.mean <- 0 + c(0.1, 0.1)%*%solve(matrix(c(1, 0.9, 0.9, 1), 2, byrow = TRUE))
z.var <- 1 - c(0.1, 0.1)%*%solve(matrix(c(1, 0.9, 0.9, 1), 2, byrow = TRUE))%*%matrix(c(0.1, 0.1), 2)

#Question 2
x <- 0
y <- 0
z <- 0
for(i in 2:1000) {
  x[i] = rnorm(1, x.mean[1]*y[i-1] + x.mean[2]*z[i-1], sqrt(x.var))
  y[i] = rnorm(1, y.mean[1]*x[i] + y.mean[2]*z[i-1], sqrt(y.var))
  z[i] = rnorm(1, z.mean[1]*x[i] + z.mean[2]*y[i], sqrt(z.var))
}

data.regular <- data.table(iteration = seq(1, 1000), x = x, method = "regular")
ggplot(data.regular, aes(x = iteration, y = x)) + geom_path() + labs(title = "Trace plot for X with regular method")

#Question 3
xy.mean <- matrix(c(0, 0), 2) + matrix(c(0.1, 0.1), 2)
xy.var <- matrix(c(1, 0.9, 0.9, 1), 2, byrow = TRUE) - matrix(c(0.1, 0.1), 2)%*%c(0.1, 0.1)

x.b <- 0
y.b <- 0
z.b <- 0
for(i in 2:1000) {
  xy = mvrnorm(1, matrix(c(0.1, 0.1), 2)*z.b[i-1], matrix(c(0.99, 0.89, 0.89, 0.99), 2, byrow = TRUE))
  x.b[i] = xy[1]
  y.b[i] = xy[2]
  z.b[i] = rnorm(1, z.mean[1]*x.b[i] + z.mean[2]*y.b[i], sqrt(z.var))
}

data.block <- data.table(iteration = seq(1, 1000), x = x.b, method = "blocked")
ggplot(data.block, aes(x = iteration, y = x)) + geom_path() + labs(title = "Trace plot for X with blocked method")

#Trace plots on the same plot
data <- rbindlist(list(data.regular, data.block))
ggplot(data, aes(x = iteration, y = x, group = method)) + geom_path(aes(color = method))
