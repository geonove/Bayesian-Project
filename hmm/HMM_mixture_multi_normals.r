# This code is combining the following two papers
# Bayesian Methods for Hidden Markov Models: Recursive Computing in the 21st Century - Author(s): Steven L. Scott
# Gibbs Sampling for Bayesian Finite Mixtures of Normal Distributions - Cristina Mollica, Luca Tardella

library("mvtnorm")

# PARAMETERS
C <- 3 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
p <- 2
mu_0 <- rep(0, p)
tau_0 <- diag(rep(0.02, p)) 
a_0 <- 5
b_0 <- diag(rep(0.5, p)) 

LENGTH <- 1000 # length of the chain
niter <- 5 # number of Gibbs Sampling iterations


# Defining the unknown mixture
w_real <- gtools::rdirichlet(1, alpha_0)
mu_real <- matrix(, nrow = C, ncol = p)
mu_real <- rmvnorm(C, mu_0, solve(tau_0))
tau_real <- array(0, dim = c(p, p, C))
tau_real[, , ] <- rWishart(C, a_0, b_0)
cat("\nw_real :", w_real, "\n")
cat("\nmu_real\n")
print(mu_real)
#cat("mu_real :", mu_real, "\n")
cat("\ntau_real\n")
print(tau_real)
#cat("tau_real :", tau_real, "\n")

# Transition probability matrix 
# For finite mixture models, all rows are equal
Q_real <- c()
for (i in 1:C) {
  Q_real <- rbind(Q_real, gtools::rdirichlet(1, alpha_0))
}
cat("\nTransition matrix\n")
print(Q_real)


# Generating an Hidden Markov Chain
HMC <- function(Q) {
  h <- c(1)
  for (i in 2:LENGTH) {
    h <- c(h, sample(1:C, prob = Q[h[length(h)], ], size = 1))
  }
  return(h)
}
h_real <- HMC(Q_real)


# Sampling from the unknown distribution
rmix <- function(w_real, mu_real, tau_real, h_real) {
  d <- matrix(, nrow = LENGTH, ncol = p)
  for (i in 1:LENGTH) {
    d[i, ] <- rmvnorm(1, mu_real[h_real[i], ], solve(tau_real[, , h_real[i]]))
  }
  
  return(d)
}
d <- rmix(w_real, mu_real, tau_real, h_real)



d <- data.frame(d)
h_real <- unlist(h_real)
x11()
plot(d, type="p", col=h_real, pch=20)



# Full conditionals
sample_Q <- function(alpha_0, h) {
  Q <- matrix(, nrow = C, ncol = C)
  NN <- matrix(0, nrow = C, ncol = C)
  for (i in 2:LENGTH) {
    NN[h[i - 1], h[i]] <- NN[h[i - 1], h[i]] + 1
  }

  for (i in 1:C) {
    Q[i, ] <- gtools::rdirichlet(1, alpha_0 + NN[i, ])
  }

  return(Q)
}


sample_tau <- function(mu, z, d, c_0, C_0, N) {
  tau <- array(, dim = c(p, p, C))

  for (c in 1:C) {
    summation <- t(sweep(d[z == c, ], p, mu[c, ], "-")) %*% (sweep(d[z == c, ], p, mu[c, ], "-"))
    c_new <- c_0 + N[c] / 2
    C_new <- C_0 + summation / 2
    tau <- rWishart(C, c_new, solve(C_new))
  }
  return(tau)
}


sample_mu <- function(tau, z, d, b_0, B_0, N) {
  mu <- matrix(, nrow = C, ncol = p)
  for (c in 1:C) {
    B_new <- solve(solve(B_0) + N[c] * tau[, , c])
    #avoid divison by 0
    if (N[c] == 0) {
      N[c] <- 1
    }
    b_new <- (B_new) %*% (solve(B_0) %*% b_0 + N[c] * (tau[, , c] %*% apply(d[z == c,], p, sum)/N[c]))
    mu[c, ] <- rmvnorm(1, b_new, B_new)
  }
  return(mu)
}


# Forward-Backward 
sample_h <- function(d, Q, mu, tau) {
  h <- numeric(LENGTH)
  # Forward recursion
  P <- array(0, dim = c(C, C, LENGTH))
  pi <- matrix(0, nrow = LENGTH, ncol = C)
  pi[1, 1] <- 1
  for (t in 2:LENGTH) {
    for (r in 1:C) {
      for (s in 1:C) {
        P[r, s, t] <- exp(log(pi[t - 1, r]) + log(Q[r, s]) + log(dmvnorm(d[t, ], mu[s,], solve(tau[, , s]))))   # dnorm -> dmvtnorm (libreria mvtfast)
      }
    }
    summation <- sum(P[, , t])

    # to reconcile proportionality
     # to avoid division by 0
    if (summation == 0) {
      P[, , t] <- matrix(1 / C^2, nrow = C, ncol = C)
    } else {
      P[, , t] <- P[, , t] / summation
    }

    for (s in 1:C) {
      pi[t, s] <-  sum(P[, s, t])
    }
  }

  # Backward recursion
  h[LENGTH] <- sample(1:C, prob = pi[LENGTH, ], size = 1)
  for (i in (LENGTH - 1):1) {
    # to avoid division by 0
    if (sum(P[, h[i + 1], i + 1]) == 0) {
      prob <- rep(1 / C, C)
    } else {
      prob <- P[, h[i + 1], i + 1]
    }
    h[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(h)
}


# Gibbs Sampler
gibbs <- function(d, niter, alpha_0, mu_0, tau_0, a_0, b_0) {
  d <- data.matrix(d)

  # Hyperparameters
  c_0 <- 2.5 * 0.5 * (p - 1)
  c_0 <- 2
  phi <- 1
  C_0 <- phi * cov(d)
  b_0 <- apply(d, p, median)
  B_0 <- diag((apply(d, p, max) - apply(d, p, min))^2)

  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) * 0.01
  Q <- c()
  for (h in 1:C) {
    Q <- rbind(Q, w)
  }
  mu <- matrix(, nrow = C, ncol = p)
  tau <- array(0, dim = c(p, p, C))
  tau[, , ]  <- rWishart(C, c_0, solve(C_0))
  mu <- rmvnorm(C, b_0, B_0)
  h <- HMC(Q)

  # Save the Markov Chain of mu
  mu_GS <- array(dim = c(C, p, niter))
  mu_GS[, , 1] <- mu

  tau_GS <- array(dim = c(p, p, C, niter))
  tau_GS[, , , 1] <- tau

  Q_GS <- array(dim = c(C, C, niter) )
  Q_GS[, , 1] <- Q
    
  h_GS <- array(dim = c(LENGTH , niter))
  h_GS[, 1] <- h

  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 1 == 0) {
      cat(i, "/", niter, "\n")
    }
    z_one_hot <- matrix(0, nrow = LENGTH, ncol = C)
    for (l in 1:LENGTH) {
      z_one_hot[l, h[l]] <- 1
    }
    N <- c()
    for (c in 1:C) {
      N <- c(N, sum(z_one_hot[, c]))
    }
    Q <- sample_Q(alpha_0, h)
    tau <- sample_tau(mu, h, d, c_0, C_0, N)
    mu <- sample_mu(tau, h, d, b_0, B_0, N)
    h <- sample_h(d, Q, mu, tau)
    mu_GS[, , i] <- mu
  }
  print(mu)
  print(Q)
  return(mu_GS)
}

mu_GS <- gibbs(d, niter, alpha_0, mu_0, tau_0, a_0, b_0)
###
#x11()
#matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)