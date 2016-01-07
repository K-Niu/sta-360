library(data.table)
library(ggplot2)
library(MASS)
library(MCMCpack)

##Reading, normalizing, and cleaning data
data <- read.csv('BP.csv', sep = ' ')
data = data.table(data)
data[, y := BP_syst/BP_dia]
data$age = (data$age - mean(data$age))/sd(data$age)
data$y = (data$y - mean(data$y))/sd(data$y)
data$id = as.factor(data$id)
ids <- levels(data$id)
data$id.new <- sapply(data$id, function(x) {
  return(which(ids == x))
})
data$id.new = as.factor(data$id.new)
n <- length(ids)
N <- length(data$id)

##Approach a
#Gibbs sampling
alpha.Gibbs <- list(lapply(seq(1, n), function(x) {
  return(matrix(c(1, 1), 2))
}))
beta.Gibbs <- list(matrix(c(1, 1), 2))
sigma.sq.inv.Gibbs <- 1
omega.Gibbs <- list(diag(2))
double.sum.xx <- Reduce('+', lapply(seq(1, n), function(i) {
  data.i = data[id.new == as.character(i),]
  x = lapply(data.i$age, function(age) {
    return(matrix(c(1, age), 2))
  })
  sum.xx = Reduce('+', lapply(x, function(xij) {
    return(xij %*% t(xij))
  }))
}))
double.sum.xy <- Reduce('+', lapply(seq(1, n), function(i) {
  data.i = data[id.new == as.character(i),]
  x = lapply(data.i$age, function(age) {
    return(matrix(c(1, age), 2))
  })
  y = data.i$y
  sum.xy = Reduce('+', lapply(seq(1, length(x)), function(j) {
    return(x[[j]]*y[j])
  }))
}))
for(k in 2:1000) {
  if(k %% 5 == 0) {
    print(k)
  }
  alpha.Gibbs[[k]] = lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    y = data.i$y
    sum.xx = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(xij))
    }))
    sum.xy = Reduce('+', lapply(seq(1, length(x)), function(j) {
      return(x[[j]]*y[j])
    }))
    sum.xbetax = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(beta.Gibbs[[k-1]]) %*% xij)
    }))
    return(mvrnorm(1, solve(sigma.sq.inv.Gibbs[k-1]*sum.xx + solve(omega.Gibbs[[k-1]])) %*% (sigma.sq.inv.Gibbs[k-1]*sum.xy - sigma.sq.inv.Gibbs[k-1]*sum.xbetax), solve(sigma.sq.inv.Gibbs[k-1]*sum.xx + solve(omega.Gibbs[[k-1]]))))
  })
  double.sum.xalphax = Reduce('+', lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    sum.xalphax = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(alpha.Gibbs[[k]][[i]]) %*% xij)
    }))
  }))
  beta.Gibbs[[k]] = mvrnorm(1, solve(sigma.sq.inv.Gibbs[k-1]*double.sum.xx + diag(2)) %*% (sigma.sq.inv.Gibbs[k-1]*double.sum.xy - sigma.sq.inv.Gibbs[k-1]*double.sum.xalphax), solve(sigma.sq.inv.Gibbs[k-1]*double.sum.xx + diag(2)))
  double.sum.SSE = Reduce('+', lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    y = data.i$y
    sum.SSE = Reduce('+', lapply(seq(1, length(x)), function(j) {
      return((y[j] - (t(x[[j]]) %*% beta.Gibbs[[k]] + t(x[[j]]) %*% alpha.Gibbs[[k]][[i]]))^2)
    }))
  }))
  sigma.sq.inv.Gibbs[k] = rgamma(1, N/2 + 1, 1 + 0.5*double.sum.SSE)
  sum.alpha = Reduce('+', lapply(alpha.Gibbs[[k]], function(alphai) {
    return(alphai %*% t(alphai))
  }))
  omega.Gibbs[[k]] = riwish(3 + n, diag(2) + sum.alpha)
}

##Approach b
#Gibbs sampling
alpha.Gibbs2 <- list(lapply(seq(1, n), function(x) {
  return(matrix(c(1, 1), 2))
}))
lambda.Gibbs2 <- list(matrix(c(1, 1), 2))
beta.Gibbs2 <- list(matrix(c(1, 1), 2))
sigma.sq.inv.Gibbs2 <- 1
for(k in 2:1000) {
  if(k %% 5 == 0) {
    print(k)
  }
  alpha.Gibbs2[[k]] = lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    y = data.i$y
    sum.xx = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(xij))
    }))
    sum.xy = Reduce('+', lapply(seq(1, length(x)), function(j) {
      return(x[[j]]*y[j])
    }))
    sum.xbetax = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(beta.Gibbs2[[k-1]]) %*% xij)
    }))
    return(mvrnorm(1, solve(diag(as.vector(lambda.Gibbs2[[k-1]])) %*% (sigma.sq.inv.Gibbs2[k-1]*sum.xx) %*% diag(as.vector(lambda.Gibbs2[[k-1]])) + diag(2)) %*% (diag(as.vector(lambda.Gibbs2[[k-1]])) %*% (sigma.sq.inv.Gibbs2[k-1]*sum.xy) - diag(as.vector(lambda.Gibbs2[[k-1]])) %*% (sigma.sq.inv.Gibbs2[k-1]*sum.xbetax)), solve(diag(as.vector(lambda.Gibbs2[[k-1]])) %*% (sigma.sq.inv.Gibbs2[k-1]*sum.xx) %*% diag(as.vector(lambda.Gibbs2[[k-1]])) + diag(2))))
  })
  lambda.Gibbs2[[k]] = lapply(seq(1, n), function(i) {
    double.sum.xx2 = Reduce('+', lapply(seq(1, n), function(i2) {
      data.i2 = data[id.new == as.character(i2),]
      x = lapply(data.i2$age, function(age) {
        return(matrix(c(1, age), 2))
      })
      sum.xx = Reduce('+', lapply(x, function(xij) {
        return(diag(as.vector(alpha.Gibbs2[[k]][[i]])) %*% xij %*% t(xij) %*% diag(as.vector(alpha.Gibbs2[[k]][[i]])))
      }))
    }))
    double.sum.xy2 = Reduce('+', lapply(seq(1, n), function(i2) {
      data.i2 = data[id.new == as.character(i2),]
      x = lapply(data.i2$age, function(age) {
        return(matrix(c(1, age), 2))
      })
      y = data.i2$y
      sum.xy = Reduce('+', lapply(seq(1, length(x)), function(j) {
        return(diag(as.vector(alpha.Gibbs2[[k]][[i]])) %*% (x[[j]]*y[j]))
      }))
    }))
    return(mvrnorm(1, solve((sigma.sq.inv.Gibbs2[k-1]*double.sum.xx2) + diag(2)) %*% ((sigma.sq.inv.Gibbs2[k-1]*double.sum.xy2) - (sigma.sq.inv.Gibbs2[k-1]*double.sum.xx2 %*% beta.Gibbs2[[k-1]])), solve((sigma.sq.inv.Gibbs2[k-1]*double.sum.xx2) + diag(2))))
  })
  double.sum.xalphax = Reduce('+', lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    sum.xalphax = Reduce('+', lapply(x, function(xij) {
      return(xij %*% t(alpha.Gibbs2[[k]][[i]]*lambda.Gibbs2[[k]]) %*% xij)
    }))
  }))
  beta.Gibbs2[[k]] = mvrnorm(1, solve(sigma.sq.inv.Gibbs2[k-1]*double.sum.xx + diag(2)) %*% (sigma.sq.inv.Gibbs2[k-1]*double.sum.xy - sigma.sq.inv.Gibbs2[k-1]*double.sum.xalphax), solve(sigma.sq.inv.Gibbs2[k-1]*double.sum.xx + diag(2)))
  double.sum.SSE = Reduce('+', lapply(seq(1, n), function(i) {
    data.i = data[id.new == as.character(i),]
    x = lapply(data.i$age, function(age) {
      return(matrix(c(1, age), 2))
    })
    y = data.i$y
    sum.SSE = Reduce('+', lapply(seq(1, length(x)), function(j) {
      return((y[j] - (t(x[[j]]) %*% beta.Gibbs2[[k]] + t(x[[j]]) %*% lambda.Gibbs2[[k]]*alpha.Gibbs2[[k]][[i]]))^2)
    }))
  }))
  sigma.sq.inv.Gibbs2[k] = rgamma(1, N/2 + 1, 1 + 0.5*double.sum.SSE)
}

##Visualizing results
beta.estimate <- Reduce('+', beta.Gibbs)/1000
visualize.group <- function(i) {
  alpha.estimate = Reduce('+', lapply(alpha.Gibbs, function(x) {
    return(x[[i]])
  }))/1000
  y.hat = sapply(data[id.new == as.character(i),]$age, function(age) {
    x = matrix(c(1, age), 2)
    return(t(x) %*% beta.estimate + t(x) %*% alpha.estimate)
  })
  estimates = data.table(age = data[id.new == as.character(i),]$age, y.hat = y.hat)
  ggplot(data[id.new == as.character(i),], aes(x = age, y = y)) + geom_point() + geom_line(data = estimates, aes(x = age, y = y.hat), color = 'royalblue1') + labs(title = paste0('Estimates for group ', i))
}
