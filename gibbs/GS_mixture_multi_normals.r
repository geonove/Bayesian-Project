# Implementing the following paper:
# https://epub.jku.at/obvulihs/download/pdf/5554146?originalFilename=true
set.seed(42)
library("mvtnorm")
# PARAMETERS
C <- 4 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
p <- 2
mu_0 <- rep(0, p)
tau_0 <- diag(rep(0.02, p)) 
a_0 <- 5
b_0 <- diag(rep(0.5, p)) 


N_SAMPLES <- 500 # number of samples
niter <- 400 # number of Gibbs Sampling iterations


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


# Sampling from the unknown distribution
rmix <- function(w_real, mu_real, tau_real) {
  z_real <- sample(1:C, prob = w_real, size = N_SAMPLES, replace = TRUE)
  x <- matrix(, nrow = N_SAMPLES, ncol = p)
  for (i in 1:N_SAMPLES) {
    x[i, ] <- rmvnorm(1, mu_real[z_real[i], ], solve(tau_real[, , z_real[i]]))
  }
  
  return(list(x, z_real))
}
mix <- rmix(w_real, mu_real, tau_real)
x <- mix[1]
z_real <- mix[2]

x <- data.frame(x)
z_real <- unlist(z_real)
x11()
plot(x, type="p", col=z_real, pch=20)


# Full conditionals
sample_w <- function(alpha_0, N) {
  alpha_new <- alpha_0 + N
  w <- gtools::rdirichlet(1, alpha_new)
  return(w)
}

sample_tau <- function(mu, z, x, c_0, C_0, N) {
  tau <- array(0, dim = c(p, p, C))

  for (c in 1:C) {
    summation <- t(sweep(x[z == c, ], p, mu[c, ], "-")) %*% (sweep(x[z == c, ], p, mu[c, ], "-"))
    c_new <- c_0 + N[c] / 2
    C_new <- C_0 + summation / 2
    tau <- rWishart(C, c_new, solve(C_new))
  }
  return(tau)
}


sample_mu <- function(tau, z, x, b_0, B_0, N) {
  mu <- matrix(0, nrow = C, ncol = p)
  for (c in 1:C) {
    B_new <- solve(solve(B_0) + N[c] * tau[, , c])
    # avoid divison by 0
    if (N[c] == 0) {
      N[c] <- 1
    }
    b_new <- (B_new) %*% (solve(B_0) %*% b_0 + N[c] * (tau[, , c] %*% apply(x[z == c,], p, sum)/N[c]))
    mu[c, ] <- rmvnorm(1, b_new, B_new)
  }
  return(mu)
}

sample_z <- function(mu, tau, w, x) {
  z <- array(dim = N_SAMPLES)
  for (i in 1:N_SAMPLES) {
    prob <- c()
    summation <- 0
    for (c in 1:C) {
      summation <- summation + w[c] * dmvnorm(x[i, ], mu[c, ], solve(tau[, , c]))
    }
    # avoid division by 0
    if (summation != 0) {
      for (c in 1:C) {
        prob <- c(prob, w[c] * dmvnorm(x[i, ], mu[c, ], solve(tau[, , c])) / summation)
      }
    } else {
      prob <- runif(n = C, min = 0, max = 1)
    }
    z[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(z)
}


# Gibbs Sampler
gibbs <- function(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0) {
  x <- data.matrix(x)

  # Hyperparameters
  c_0 <- 2.5 * 0.5 * (p - 1)
  c_0 <- 2
  phi <- 1
  C_0 <- phi * cov(x)
  b_0 <- apply(x, p, median)
  B_0 <- diag((apply(x, p, max) - apply(x, p, min))^2)
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) * 0.01

  mu <- matrix(, nrow = C, ncol = p)
  tau <- array(0, dim = c(p, p, C))
  tau[, , ]  <- rWishart(C, c_0, solve(C_0))
  mu <- rmvnorm(C, b_0, B_0)
  z <- array(dim = N_SAMPLES)
  for (i in 1:N_SAMPLES) {
    z[i] <- sample(1:C, prob = w, size = 1)
  }

  # Save the Markov Chain of mu
  mu_GS <- array(dim = c(C, p, niter))
  mu_GS[, , 1] <- mu

  tau_GS <- array(dim = c(p, p, C, niter))
  tau_GS[, , , 1] <- tau

  w_GS <- array(dim = c(C, niter) )
  w_GS[, 1] <- w


    
  z_GS <- array(dim = c(N_SAMPLES , niter))
  z_GS[, 1] <- z

  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 10 == 0) {
      cat(i, "/", niter, "\n")
    }

    N <- as.data.frame(table(z))$Freq

    w <- sample_w(alpha_0, N)
    tau <- sample_tau(mu, z, x, c_0, C_0, N)
    mu <- sample_mu(tau, z, x, b_0, B_0, N)
    z <- sample_z(mu, tau, w, x)

    mu_GS[, , i] <- mu
    tau_GS[, , , i] <- tau
    w_GS[, i] <- w
    z_GS[, i] <- z
  }

  return(list("mu_GS" = mu_GS, "tau_GS" = tau_GS, "w_GS" = w_GS, "z_GS" = z_GS))
}

mix <- gibbs(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0)
mu_GS <- mix$mu_GS

for ( i in 1:p ) {

  x11()
  matplot(t(mu_GS[, i, ]), main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)
}
z_GS <- data.frame(mix[4])
x11()
plot(x, type="p", col=z_GS[, niter], pch=20)
