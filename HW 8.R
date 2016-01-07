library(MASS)
library(glmnet)
library(data.table)
library(ggplot2)

##Generate data
n = 100
beta.true <- c(rnorm(5, 0, sqrt(2)), rep(0, 95))
beta <- matrix(beta.true, n)
x <- mvrnorm(n, matrix(rep(0, n), n), diag(n))
y <- beta.true[1]*x[, 1] + beta.true[2]*x[, 2] + beta.true[3]*x[, 3] + beta.true[4]*x[, 4] + beta.true[5]*x[, 5] + rnorm(n, 0, 1)

##Fully Bayes ridge
#Prior
a <- 0.5
b <- 0.5
c <- 0.001
d <- 0.001
#Gibbs sampling
beta.ridge <- list(beta)
tau.ridge <- 1
sigma.sq.inv.ridge <- 1
for(i in 2:1000) {
  beta.ridge[[i]] = matrix(mvrnorm(1, solve(t(x) %*% x + tau.ridge[i-1]*diag(100)) %*% (t(x) %*% y), solve(t(x) %*% x + tau.ridge[i-1]*diag(100))*(1/sigma.sq.inv.ridge[i-1])), 100)
  tau.ridge[i] = rgamma(1, a + 50, rate = b + 0.5*sigma.sq.inv.ridge[i-1]*(t(beta.ridge[[i]]) %*% beta.ridge[[i]]))
  sigma.sq.inv.ridge[i] = rgamma(1, c + 100, d + 0.5*(t(y - x %*% beta.ridge[[i]]) %*% (y - x %*% beta.ridge[[i]])) + 0.5*tau.ridge[i]*(t(beta.ridge[[i]]) %*% beta.ridge[[i]]))
}
beta.ridge.matrix <- matrix(unlist(beta.ridge), ncol = 100, byrow = TRUE)
ridge.beta <- colMeans(beta.ridge.matrix)

##t-prior with v = 1
#Prior
a <- 0.5
b <- 0.5
c <- 0.001
d <- 0.001
v <- 1
#Gibbs sampling
beta.t <- list(beta)
lambda.t <- list(rep(1, 100))
tau.t <- 1
sigma.sq.inv.t <- 1
for(i in 2:1000) {
  lambda.t[[i]] = sapply(seq(1, 100), function(j, v, beta, tau) {
    return(rgamma(1, 0.5*(v + 1), rate = 0.5*(tau*beta[j]^2 + v)))
  }, v = v, beta = beta.t[[i-1]], tau = tau.t[i-1])
  lambda.matrix = diag(lambda.t[[i]])
  beta.t[[i]] = matrix(mvrnorm(1, solve(sigma.sq.inv.t[i-1]*(t(x) %*% x) + tau.t[i-1]*lambda.matrix) %*% (sigma.sq.inv.t[i-1]*(t(x) %*% y)), solve(sigma.sq.inv.t[i-1]*(t(x) %*% x) + tau.t[i-1]*lambda.matrix)), 100)
  tau.t[i] = rgamma(1, a + 50, rate = b + 0.5*(t(beta.t[[i]]) %*% lambda.matrix %*% beta.t[[i]]))
  sigma.sq.inv.t[i] = rgamma(1, c + 50, d + 0.5*(t(y - x %*% beta.t[[i]]) %*% (y - x %*% beta.t[[i]])))
}
beta.t.1.matrix <- matrix(unlist(beta.t), ncol = 100, byrow = TRUE)
t.1.beta <- colMeans(beta.t.1.matrix)

##t-prior with v = 0.001
#Prior
a <- 0.5
b <- 0.5
c <- 0.001
d <- 0.001
v <- 0.001
#Gibbs sampling
beta.t <- list(beta)
lambda.t <- list(rep(1, 100))
tau.t <- 1
sigma.sq.inv.t <- 1
for(i in 2:1000) {
  lambda.t[[i]] = sapply(seq(1, 100), function(j, v, beta, tau) {
    return(rgamma(1, 0.5*(v + 1), rate = 0.5*(tau*beta[j]^2 + v)))
  }, v = v, beta = beta.t[[i-1]], tau = tau.t[i-1])
  lambda.matrix = diag(lambda.t[[i]])
  beta.t[[i]] = matrix(mvrnorm(1, solve(sigma.sq.inv.t[i-1]*(t(x) %*% x) + tau.t[i-1]*lambda.matrix) %*% (sigma.sq.inv.t[i-1]*(t(x) %*% y)), solve(sigma.sq.inv.t[i-1]*(t(x) %*% x) + tau.t[i-1]*lambda.matrix)), 100)
  tau.t[i] = rgamma(1, a + 50, rate = b + 0.5*(t(beta.t[[i]]) %*% lambda.matrix %*% beta.t[[i]]))
  sigma.sq.inv.t[i] = rgamma(1, c + 50, d + 0.5*(t(y - x %*% beta.t[[i]]) %*% (y - x %*% beta.t[[i]])))
}
beta.t.001.matrix <- matrix(unlist(beta.t), ncol = 100, byrow = TRUE)
t.001.beta <- colMeans(beta.t.001.matrix)

##Lasso
model.lasso.cv <- cv.glmnet(x, y, alpha = 1)
index <- which(model.lasso.cv$lambda == model.lasso.cv$lambda.min)
lasso.beta <- model.lasso.cv$glmnet.fit$beta[,index]

##Comparison between methods
data.true <- data.table(num = seq(1, 100), beta = beta.true, group = 'true')
data.ridge <- data.table(num = seq(1, 100), beta = ridge.beta, group = 'ridge')
data.t.1 <- data.table(num = seq(1, 100), beta = t.1.beta, group = 't, v = 1')
data.t.001 <- data.table(num = seq(1, 100), beta = t.001.beta, group = 't, v = 0.001')
data.lasso <- data.table(num = seq(1, 100), beta = lasso.beta, group = 'Lasso')
data <- rbindlist(list(data.true, data.ridge, data.t.1, data.t.001, data.lasso))
ggplot(data, aes(x = num, y = beta, fill = group)) + geom_point(aes(color = group, shape = group)) + labs(title = 'Beta estimates from different methods') + scale_shape(solid = FALSE)

#MSE on betas
mean((beta.true - ridge.beta)^2)
mean((beta.true - t.1.beta)^2)
mean((beta.true - t.001.beta)^2)
mean((beta.true - lasso.beta)^2)

#Predictive MSE
x.test <- mvrnorm(100, matrix(rep(0, 100), 100), diag(100))
y.test <- beta.true[1]*x.test[, 1] + beta.true[2]*x.test[, 2] + beta.true[3]*x.test[, 3] + beta.true[4]*x.test[, 4] + beta.true[5]*x.test[, 5] + rnorm(100, 0, 1)
mean((y.test - x.test %*% ridge.beta)^2)
mean((y.test - x.test %*% t.1.beta)^2)
mean((y.test - x.test %*% t.001.beta)^2)
mean((y.test - x.test %*% lasso.beta)^2)
