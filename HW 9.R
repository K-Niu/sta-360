##Setup
n.groups <- 10
n.people <- round(runif(n.groups, 1, 20))
t <- lapply(n.people, function(n.i) {
  return(round(runif(n.i, 1, 120)))
})
phi <- rgamma(1, 2, rate = 3)
beta <- rgamma(n.groups, phi, rate = phi)
lambda <- rgamma(1, 4, rate = 5)
y <- lapply(seq(1, n.groups), function(i, t, beta, lambda, n.people) {
  return(sapply(seq(1, n.people[i]), function(j, i, t, beta, lambda) {
    return(rpois(1, t[[i]][[j]]*beta[i]*lambda))
  }, i = i, t = t, beta = beta, lambda = lambda))
}, t = t, beta = beta, lambda = lambda, n.people = n.people)

##Sampling from joint posterior
#Priors
a <- 0.001
b <- 0.001
c <- 0.001
d <- 0.001
#Combination of Gibbs and Metropolis with normal random walk
beta.Gibbs <- list(rep(1, n.groups))
lambda.Gibbs <- 1
phi.Metro <- 1
likelihood <- function(beta, lambda, y, t, n.groups, n.people) {
  return(sum(sapply(seq(1, n.groups), function(i, beta, lambda, y, t, n.people) {
    return(sum(sapply(seq(1, n.people[i]), function(j, i, beta, lambda, y, t) {
      return(dpois(y[[i]][[j]], t[[i]][[j]]*beta[i]*lambda, log = TRUE))
    }, i = i, beta = beta, lambda = lambda, y = y, t = t)))
  }, beta = beta, lambda = lambda, y = y, t = t, n.people = n.people)))
}
priors <- function(beta, lambda, phi, a, b, c, d, n.groups) {
  beta.prior = sum(sapply(seq(1, n.groups), function(i, beta, phi) {
    return(dgamma(beta[i], phi, rate = phi, log = TRUE))
  }, beta = beta, phi = phi))
  lambda.prior = dgamma(lambda, a, rate = b, log = TRUE)
  phi.prior = dgamma(phi, c, rate = d, log = TRUE)
  return(beta.prior + lambda.prior + phi.prior)
}
f <- function(beta, lambda, phi, y, t, a, b, c, d, n.groups, n.people) {
  return(likelihood(beta, lambda, y, t, n.groups, n.people) + priors(beta, lambda, phi, a, b, c, d, n.groups))
}
accept.count <- 0
for(k in 2:10000) {
  beta.Gibbs[[k]] = sapply(seq(1, n.groups), function(i, lambda, phi, y, t) {
    rgamma(1, phi + sum(y[[i]]), rate = phi + lambda*sum(t[[i]]))
  }, lambda = lambda.Gibbs[k-1], phi = phi.Metro[k-1], y = y, t = t)
  lambda.Gibbs[[k]] = rgamma(1, a + sum(sapply(y, sum)), rate = b + sum(sapply(seq(1, n.groups), function(i, beta, t) {
    return(beta[i]*sum(t[[i]]))
  }, beta = beta.Gibbs[[k]], t = t)))
  #Metropolis
  phi.cand = rnorm(1, phi.Metro[k-1], sqrt(0.1))
  f.cand = f(beta.Gibbs[[k]], lambda.Gibbs[k], phi.cand, y, t, a, b, c, d, n.groups, n.people)
  f.cur = f(beta.Gibbs[[k]], lambda.Gibbs[k], phi.Metro[k-1], y, t, a, b, c, d, n.groups, n.people)
  if(is.na(f.cand)) {
    phi.Metro[k] = phi.Metro[k-1]
    #print('Reject')
  } else {
    prob = min(1, exp(f.cand - f.cur))
    accept = rbinom(1, 1, prob)
    if(accept == 1) {
      phi.Metro[k] = phi.cand
      #print('Accept')
      accept.count = accept.count + 1
    } else {
      phi.Metro[k] = phi.Metro[k-1]
      #print('Reject')
    }
  }
}
#Trace plots
beta.matrix <- matrix(unlist(beta.Gibbs), ncol = n.groups, byrow = TRUE)
plot(seq(1, 10000), beta.matrix[,1], type = 'l', main = 'Trace plot for beta1', xlab = 'iteration', ylab = 'beta1')
plot(seq(1, 10000), lambda.Gibbs, type = 'l', main = 'Trace plot for lambda', xlab = 'iteration', ylab = 'lambda')
plot(seq(1, 10000), phi.Metro, type = 'l', main = 'Trace plot for phi', xlab = 'iteration', ylab = 'phi')
