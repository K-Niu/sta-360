library(data.table)
library(ggplot2)
library(MASS)
library(MCMCpack)

data <- fread('HW 11 - data.csv')
data[, Y := (Y - mean(Y))/sqrt(var(Y))]
data[, missing := ifelse(!is.na(X.1) & !is.na(X.2), FALSE, TRUE)]
data.obs <- data[missing == FALSE,]
data.mis <- data[missing == TRUE,] 
yobs <- data.obs$Y
xobs <- as.matrix(data.obs)[,2:3]
xobs = t(xobs)
ymis <- data.mis$Y
xmis <- as.matrix(data.mis)[,2:3]
xmis = t(xmis)
n <- length(yobs) + length(ymis)

#Priors
m.b <- matrix(c(1, 1), 2)
s.b <- diag(2)
a <- 0.001
b <- 0.001
m.v <- matrix(c(1, 1), 2)
s.v <- diag(2)
#Gibbs
xmis.init <- xmis
xmis.init[is.na(xmis.init)] = 1
xmis.Gibbs <- list(xmis.init)
beta.Gibbs <- list(matrix(c(1, 1), 2))
mu.Gibbs <- list(matrix(c(1, 1), 2))
sigma.Gibbs <- list(diag(2))
tau.Gibbs <- 1
sum.xx.obs <- Reduce('+', lapply(seq(1, ncol(xobs)), function(i) {
  return(xobs[,i] %*% t(xobs[,i]))
}))
sum.xy.obs <- Reduce('+', lapply(seq(1, ncol(xobs)), function(i) {
  return(xobs[,i]*yobs[i])
}))
sum.x.obs <- Reduce('+', lapply(seq(1, ncol(xobs)), function(i) {
  return(xobs[,i])
}))
for(k in 2:10000) {
  xmis.temp = xmis
  xmis.Gibbs[[k]] = sapply(seq(1, ncol(xmis.temp)), function(i) {
    xmis.mu = solve(tau.Gibbs[k-1]*(beta.Gibbs[[k-1]] %*% t(beta.Gibbs[[k-1]])) + solve(sigma.Gibbs[[k-1]])) %*% (tau.Gibbs[k-1]*beta.Gibbs[[k-1]]*ymis[i] + solve(sigma.Gibbs[[k-1]]) %*% mu.Gibbs[[k-1]])
    xmis.sigma = solve(tau.Gibbs[k-1]*(beta.Gibbs[[k-1]] %*% t(beta.Gibbs[[k-1]])) + solve(sigma.Gibbs[[k-1]]))
    xmis.temp.i = xmis.temp[,i]
    if(is.na(xmis.temp.i[1]) & is.na(xmis.temp.i[2])) {
      return(mvrnorm(1, xmis.mu, xmis.sigma))
    } else if(is.na(xmis.temp.i[1])) {
      return(c(rnorm(1, xmis.mu[1] + xmis.sigma[1, 2]*(1/xmis.sigma[2, 2])*(xmis.temp.i[2] - xmis.mu[2]), sqrt(xmis.sigma[1, 1] - xmis.sigma[1, 2]*(1/xmis.sigma[2, 2])*xmis.sigma[2, 1])), xmis.temp.i[2]))
    } else {
      return(c(xmis.temp.i[1], rnorm(1, xmis.mu[2] + xmis.sigma[2, 1]*(1/xmis.sigma[1, 1])*(xmis.temp.i[1] - xmis.mu[1]), sqrt(xmis.sigma[2, 2] - xmis.sigma[2, 1]*(1/xmis.sigma[1, 1])*xmis.sigma[1, 2]))))
    }
  })
  sum.xx.mis = Reduce('+', lapply(seq(1, ncol(xmis.Gibbs[[k]])), function(i) {
    return(xmis.Gibbs[[k]][,i] %*% t(xmis.Gibbs[[k]][,i]))
  }))
  sum.xy.mis = Reduce('+', lapply(seq(1, ncol(xmis.Gibbs[[k]])), function(i) {
    return(xmis.Gibbs[[k]][,i]*ymis[i])
  }))
  sum.x.mis = Reduce('+', lapply(seq(1, ncol(xmis.Gibbs[[k]])), function(i) {
    return(xmis.Gibbs[[k]][,i])
  }))
  beta.Gibbs[[k]] = mvrnorm(1, solve(tau.Gibbs[k-1]*sum.xx.obs + tau.Gibbs[k-1]*sum.xx.mis + solve(s.b)) %*% (tau.Gibbs[k-1]*sum.xy.obs + tau.Gibbs[k-1]*sum.xy.mis + solve(s.b) %*% m.b), solve(tau.Gibbs[k-1]*sum.xx.obs + tau.Gibbs[k-1]*sum.xx.mis + solve(s.b)))
  mu.Gibbs[[k]] = mvrnorm(1, solve(length(yobs)*solve(sigma.Gibbs[[k-1]]) + length(ymis)*solve(sigma.Gibbs[[k-1]]) + solve(s.v)) %*% (solve(sigma.Gibbs[[k-1]]) %*% sum.x.obs + solve(sigma.Gibbs[[k-1]]) %*% sum.x.mis + solve(s.v) %*% m.v), solve(length(yobs)*solve(sigma.Gibbs[[k-1]]) + length(ymis)*solve(sigma.Gibbs[[k-1]]) + solve(s.v)))
  sum.xxm.obs = Reduce('+', lapply(seq(1, ncol(xobs)), function(i) {
    return((xobs[,i] - mu.Gibbs[[k]]) %*% t(xobs[,i] - mu.Gibbs[[k]]))
  }))  
  sum.xxm.mis = Reduce('+', lapply(seq(1, ncol(xmis.Gibbs[[k]])), function(i) {
    return((xmis.Gibbs[[k]][,i] - mu.Gibbs[[k]]) %*% t(xmis.Gibbs[[k]][,i] - mu.Gibbs[[k]]))
  }))
  sigma.Gibbs[[k]] = riwish(3 + n, sum.xxm.obs + sum.xxm.mis + diag(2))
  tau.Gibbs[k] = rgamma(1, a + n/2, b + 0.5*(sum((yobs - t(xobs) %*% beta.Gibbs[[k]])^2)) + 0.5*(sum((ymis - t(xmis.Gibbs[[k]]) %*% beta.Gibbs[[k]])^2)))
}
#Trace plots
ggplot(data.table(iteration = seq(1, 10000), tau = tau.Gibbs), aes(x = iteration, y = tau)) + geom_line() + labs(title = 'Trace plot for tau')
beta.1 <- sapply(beta.Gibbs, function(x) {
  return(x[1])
})
ggplot(data.table(iteration = seq(1, 10000), beta.1 = beta.1), aes(x = iteration, y = beta.1)) + geom_line() + labs(title = 'Trace plot for beta.1')
mu.1 <- sapply(mu.Gibbs, function(x) {
  return(x[1])
})
ggplot(data.table(iteration = seq(1, 10000), mu.1 = mu.1), aes(x = iteration, y = mu.1)) + geom_line() + labs(title = 'Trace plot for mu.1')
sigma.11 <- sapply(sigma.Gibbs, function(x) {
  return(x[1,1])
})
ggplot(data.table(iteration = seq(1, 10000), sigma.11 = sigma.11), aes(x = iteration, y = sigma.11)) + geom_line() + labs(title = 'Trace plot for sigma.11')
xmis1.2 <- sapply(xmis.Gibbs, function(x) {
  return(x[2,1])
})
ggplot(data.table(iteration = seq(1, 10000), xmis1.2 = xmis1.2), aes(x = iteration, y = xmis1.2)) + geom_line() + labs(title = 'Trace plot for xmis1.2')

#Visualize estimates
xmis.est <- Reduce('+', xmis.Gibbs)/10000
data.mis.est <- data.table(Y = ymis, X.1 = xmis.est[1,], X.2 = xmis.est[2,], missing = TRUE)
data.comp <- rbindlist(list(data.obs, data.mis.est))
ggplot(data.comp, aes(x = X.1, y = X.2, group = missing)) + geom_point(aes(color = missing)) + labs(title = 'Comparing nonmissing x values and estimated missing x values')
ggplot(data.comp, aes(x = Y, group = missing)) + geom_bar(aes(fill = missing), position = 'dodge') + labs(title = 'Comparing y values for nonmissing and missing')

#Comparison to MLE
beta.est <- Reduce('+', beta.Gibbs)/10000
lm.model <- glm(Y ~ X.1 + X.2 - 1, data = data)
lm.model$coefficients

bayes <- t(beta.est) %*% xobs
freq <- matrix(lm.model$coefficients, 1) %*% xobs
estimates <- data.table(i = seq(1, 50), y = data.obs$Y, y.bayes = as.vector(bayes), y.freq = as.vector(freq))
ggplot(estimates, aes(x = i, y = y, color = 'Actual')) + geom_point() + 
  geom_point(aes(x = i, y = y.bayes, color = 'Bayes')) + 
  geom_point(aes(x = i, y = y.freq, color = 'MLE')) + 
  scale_color_manual(breaks = c('Actual', 'Bayes', 'MLE'), values = c('black', 'red', 'blue')) +
  labs(title = 'Predicted y values for Bayesian and MLE methods', color = '')
bayes.MSE <- mean((yobs - bayes)^2)
freq.MSE <- mean((yobs - freq)^2)
