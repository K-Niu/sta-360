library(MASS)
library(data.table)
library(ggplot2)

#Simulates data with covariance p
#Returns a list with the first element being the x values and the second element being the y values
simulate_data <- function(p) {
  beta = matrix(c(-1, 1, 1), 3)
  err.var = 1
  n = 100
  x = mvrnorm(n, matrix(c(0, 0), 2), matrix(c(1, p, p, 1), 2, byrow = TRUE))
  y = beta[1] + beta[2]*x[,1] + beta[3]*x[,2] + rnorm(n, 0, sqrt(err.var))
  return(list(x, y))
}

#Calculates the beta hat values and MSE of the betas
calculate_MLE <- function(x, y) {
  n = 100
  x_tilda = cbind(rep(1, n), x)
  beta_hat = solve(t(x_tilda)%*%x_tilda)%*%t(x_tilda)%*%y
  MSE = mean((c(-1, 1, 1) - beta_hat)^2)
  return(c(beta_hat, MSE))
}

#Comparing beta hat values and MSE for different values of p
p <- seq(-0.99, 0.99, 0.001)
result <- sapply(p, function(x) {
  do <- function(x) {
    data = simulate_data(x)
    values = calculate_MLE(data[[1]], data[[2]])
    return(values)
  }
  return(colMeans(t(sapply(rep(x, 30), do))))
})
result = t(result)
result.dt <- data.table(p = p, beta0 = result[,1], beta1 = result[,2], beta2 = result[,3], MSE = result[,4])
ggplot(result.dt, aes(x = p, y = MSE)) + geom_point() + labs(title = "MLE MSE on betas for different values of p")

##Ridge regression with hyperprior
data <- simulate_data(0.1)
x <- data[[1]]
y <- data[[2]]
n <- 100
#Priors
a <- 0.5
b <- 0.5
c <- 0.001
d <- 0.001
#Gibbs sampling, no burn-ins dropped because starting values of parameters were close to actually value and many iterations were used
beta_0 <- -1
beta_1 <- 1
beta_2 <- 1
tau <- 1
var_rec <- 1 #Reciprocal of var
for(i in 2:10000) {
  beta_0_mean = (1/(n + tau[i-1]))*(sum(y) - beta_1[i-1]*sum(x[,1]) - beta_2[i-1]*sum(x[,2]))
  beta_0_var = 1/((n + tau[i-1])*var_rec[i-1])
  beta_0[i] = rnorm(1, beta_0_mean, sqrt(beta_0_var))
  beta_1_mean = (1/(sum(x[,1]^2) + tau[i-1]))*(sum(x[,1]*y) - beta_0[i]*sum(x[,1]) - beta_2[i-1]*sum(x[,1]*x[,2]))
  beta_1_var = 1/((sum(x[,1]^2) + tau[i-1])*var_rec[i-1])
  beta_1[i] = rnorm(1, beta_1_mean, sqrt(beta_1_var))
  beta_2_mean = (1/(sum(x[,2]^2) + tau[i-1]))*(sum(x[,2]*y) - beta_0[i]*sum(x[,2]) - beta_1[i]*sum(x[,1]*x[,2]))
  beta_2_var = 1/((sum(x[,2]^2) + tau[i-1])*var_rec[i-1])
  beta_2[i] = rnorm(1, beta_2_mean, sqrt(beta_2_var))
  tau_a = a + 3/2
  tau_b = b + (1/2)*var_rec[i-1]*(beta_0[i]^2 + beta_1[i]^2 + beta_2[i]^2)
  tau[i] = rgamma(1, tau_a, rate = tau_b)
  var_rec_a = c + (n + 3)/2
  var_rec_b = d + (1/2)*sum((y - (beta_0[i] + beta_1[i]*x[,1] + beta_2[i]*x[,2]))^2) + (tau[i]/2)*(beta_0[i]^2 + beta_1[i]^2 + beta_2[i]^2)
  var_rec[i] = rgamma(1, var_rec_a, rate = var_rec_b)
}
ggplot(data.table(tau = tau), aes(x = tau)) + geom_density() + labs(title = 'Posterior distribution of tau from Gibbs sampling')
#MSE
mean((c(-1, 1, 1) - c(mean(beta_0), mean(beta_1), mean(beta_2)))^2)

#Calculate MSE for penalty tau on data x (matrix), y
calculate_Bayes_MSE <- function(x, y, tau) {
  n = 100
  #Split into training an test set
  x.train = x[1:90,]
  y.train = y[1:90]
  x.test = x[91:100,]
  y.test = y[91:100]
  #Prior
  c = 0.001
  d = 0.001
  #Gibbs sampling
  beta_0 = -1
  beta_1 = 1
  beta_2 = 1
  var_rec = 1 #Reciprocal of var
  for(i in 2:10000) {
    beta_0_mean = (1/(n + tau))*(sum(y.train) - beta_1[i-1]*sum(x.train[,1]) - beta_2[i-1]*sum(x.train[,2]))
    beta_0_var = 1/((n + tau)*var_rec[i-1])
    beta_0[i] = rnorm(1, beta_0_mean, sqrt(beta_0_var))
    beta_1_mean = (1/(sum(x.train[,1]^2) + tau))*(sum(x.train[,1]*y.train) - beta_0[i]*sum(x.train[,1]) - beta_2[i-1]*sum(x.train[,1]*x.train[,2]))
    beta_1_var = 1/((sum(x.train[,1]^2) + tau)*var_rec[i-1])
    beta_1[i] = rnorm(1, beta_1_mean, sqrt(beta_1_var))
    beta_2_mean = (1/(sum(x.train[,2]^2) + tau))*(sum(x.train[,2]*y.train) - beta_0[i]*sum(x.train[,2]) - beta_1[i]*sum(x.train[,1]*x.train[,2]))
    beta_2_var = 1/((sum(x.train[,2]^2) + tau)*var_rec[i-1])
    beta_2[i] = rnorm(1, beta_2_mean, sqrt(beta_2_var))
    var_rec_a = c + (n + 3)/2
    var_rec_b = d + (1/2)*sum((y.train - (beta_0[i] + beta_1[i]*x.train[,1] + beta_2[i]*x.train[,2]))^2) + (tau/2)*(beta_0[i]^2 + beta_1[i]^2 + beta_2[i]^2)
    var_rec[i] = rgamma(1, var_rec_a, rate = var_rec_b)
  }
  y_hat = mean(beta_0) + mean(beta_1)*x.test[,1] + mean(beta_2)*x.test[,2]
  MSE = mean((y.test - y_hat)^2)
  return(MSE)
}

#Ridge regression with cross-validation
#Takes a long time
tau.selected <- NULL
for(i in 1:100) {
  data = simulate_data(0.1)
  x = data[[1]]
  y = data[[2]]
  tau.grid = seq(0.1, 5, 0.1)
  tau.MSE = sapply(tau.grid, calculate_Bayes_MSE, x = x, y = y)
  tau.selected[i] = tau.grid[which.min(tau.MSE)]
}
ggplot(data.table(tau = tau.selected), aes(x = tau)) + geom_histogram() + labs(title = 'Selection of tau values by cross-validation')
