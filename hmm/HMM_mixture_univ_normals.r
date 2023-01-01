# This code is combining the following two papers
# Bayesian Methods for Hidden Markov Models: Recursive Computing in the 21st Century - Author(s): Steven L. Scott
# Gibbs Sampling for Bayesian Finite Mixtures of Normal Distributions - Cristina Mollica, Luca Tardella


# PARAMETERS
C <- 3 # number of possible hidden states
alpha_0 <- runif(n = C, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

LENGTH <- 500 # length of the chain
niter <- 100 # number of Gibbs Sampling iterations


# Defining the unknown mixture
mu_real <- rnorm(C, mu_0, sqrt(1 / tau_0))
tau_real <- rgamma(C, shape = a_0, rate = b_0)
cat("\nmu_real\n")
print(mu_real)
cat("\ntau_real\n")
print(tau_real)

# Transition probability matrix 
# (for finite mixture models all rows are equal)
Q_real <- c()
for (i in 1:C) {
  Q_real <- rbind(Q_real, gtools::rdirichlet(1, alpha_0))
}
cat("\nQ_real\n")
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
rmix <- function(h_real, mu_real, tau_real) {
  d <- rnorm(length(h_real), mu_real[h_real], sqrt(1 / tau_real[h_real]))
  return(d)
}
d <- rmix(h_real, mu_real, tau_real)


# Plotting the Hidden Markov Model and samples
library(rgl)
plot_HMM_samples <- function(x_seq, d, h_real) {
  norms <- data.frame(x_seq)
  norms <- c()
  for (i in 1:C) {
    norms <- rbind(norms, dnorm(x_seq, mu_real[i], sqrt(1 / tau_real[i])))
  }

  for (i in 1:LENGTH) {
    plot3d(x_seq, norms[h_real[i], ], rep(i, length(x_seq)), type = "l", lwd = 4, col = h_real[i], zlim = c(1, LENGTH))
    plot3d(d[i], 0, i, size = 4, col = h_real[i], zlim = c(1, LENGTH))
  }
  #grid3d(c("x", "y+", "z"))
}
x_seq <- seq(min(d) - 2*max(sqrt(1/tau_real)), max(d)+ 2*max(sqrt(1/tau_real)), by = 0.001)
#plot_HMM_samples(x_seq, d, h_real)


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


sample_tau <- function(mu, z, x, a_0, b_0, N) {
  tau <- c()
  for (c in 1:C) {
    summation <- 0
    for (i in 1:LENGTH) {
      summation <- summation + z[i, c] * (x[i] - mu[c])^2
    }
    tau <- c(tau, rgamma(1, shape = a_0 + N[c] / 2, rate = b_0 + summation / 2))
  }
  return(tau)
}


sample_mu <- function(tau, z, x, tau_0, mu_0, N) {
  mu <- c()
  for (c in 1:C) {

    # to avoid division by 0
    if (N[c] == 0) {
      N[c] <- 1
    }

    delta <- N[c] * tau[c] / (tau_0 + N[c] * tau[c])
    x_bar <- (z[, c] %*% x) / N[c]
    mu <- c(mu, rnorm(1, delta * x_bar + (1 - delta) * mu_0, sqrt(1 / (tau_0 + N[c] * tau[c]))))
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
        P[r, s, t] <- exp(log(pi[t - 1, r]) + log(Q[r, s]) + log(dnorm(d[t], mu[s], sqrt(1 / tau[s]))))   # dnorm -> dmvtnorm (libreria mvtfast)
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
  h <- sample(1:C, prob = pi[LENGTH, ], size = 1)
  for (i in (LENGTH - 1):1) {
    # to avoid division by 0
    if (sum(P[, h[1], i + 1]) == 0) {
      prob <- rep(1 / C, C)
    } else {
      prob <- P[, h[1], i + 1]
    }
    h[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(h)
}


# Gibbs Sampler
gibbs <- function(d, niter, alpha_0, mu_0, tau_0, a_0, b_0) {
  b_0 <- 50

  cat("\n\nGibbs Sampler\n")
  w <- gtools::rdirichlet(1, alpha_0)
  Q <- c()
  for (h in 1:C) {
    Q <- rbind(Q, w)
  }
  tau <- rgamma(C, shape = a_0, rate = b_0)
  mu <- rnorm(C, mu_0, sqrt(1 / tau_0))
  h <- HMC(Q)

  mu_GS <- matrix(, nrow = niter, ncol = C)
  mu_GS[1, ] <- mu

  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 50 == 0) {
      cat(i, "/", niter, "\n")
    }
    z <- matrix(0, nrow = LENGTH, ncol = C)
    for (l in 1:LENGTH) {
      z[l, h[l]] <- 1
    }
    N <- c()
    for (c in 1:C) {
      N <- c(N, sum(z[, c]))
    }
    Q <- sample_Q(alpha_0, h)
    tau <- sample_tau(mu, z, d, a_0, b_0, N)
    mu <- sample_mu(tau, z, d, tau_0, mu_0, N)
    h <- sample_h(d, Q, mu, tau)
    mu_GS[i, ] <- mu
  }
  cat("\nmu\n")
  print(mu)
  cat("\nQ\n")
  print(Q)
  return(mu_GS)
}

mu_GS <- gibbs(d, niter, alpha_0, mu_0, tau_0, a_0, b_0)

x11()
matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)