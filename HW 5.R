#Set values
a <- 4
b <- 1
c <- 1
d <- 1
y <- 3
x <- 6

#Perfect sampling
theta_perfect <- NULL
gamma_perfect <- NULL
for(i in 1:1000) {
  gamma_perfect[i] = rgamma(1, c + x, d)
  theta_perfect[i] = rgamma(1, a + x + y, b + 1 + gamma_perfect[i])
}

#Gibbs sampling
theta_gibbs <- 4
gamma_gibbs <- 1
for(i in 2:1000) {
  gamma_gibbs[i] = rgamma(1, c + x, d + theta_gibbs[i-1])
  theta_gibbs[i] = rgamma(1, a + x + y, rate = b + 1 + gamma_gibbs[i])
}

#Posterior summaries
mean(theta_perfect)
mean(gamma_perfect)
mean(theta_gibbs)
mean(gamma_gibbs)
mean(gamma_gibbs > 1) #P(gamma > 1 | y, x)
theta_perfect = sort(theta_perfect)
gamma_perfect = sort(gamma_perfect)
theta_gibbs = sort(theta_gibbs)
gamma_gibbs = sort(gamma_gibbs)
c(theta_perfect[25], theta_perfect[975])
c(gamma_perfect[25], gamma_perfect[975])
c(theta_gibbs[25], theta_gibbs[975])
c(gamma_gibbs[25], gamma_gibbs[975])

#Repeat for three groups
a <- 4
b <- 1
phi <- 1
y <- 3
x <- 6
z <- 2

theta <- 4
gamma <- 1
psi <- 1
for(i in 2:1000) {
  gamma[i] = rgamma(1, phi + x, phi + theta[i-1])
  psi[i] = rgamma(1, phi + z, phi + theta[i-1])
  theta[i] = rgamma(1, a + x + y + z, b + 1 + gamma[i] + psi[i])
}

mean(theta)
mean(gamma)
mean(psi)
mean(gamma > 1)
mean(psi < 1)
mean(gamma > psi)
theta = sort(theta)
gamma = sort(gamma)
psi = sort(psi)
c(theta[25], theta[975])
c(gamma[25], gamma[975])
c(psi[25], psi[975])