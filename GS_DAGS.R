library(gRbase)
library(GLSME)
library(BCDAG)
library(gss)
library(graph)
library(Rcpp)
library(rgl)
library(car)

setwd("C:/Users/andre/Documents/Università/Magistrale/2° Anno/1° Semestre/Bayesian statistics/Progetto/codes")
#setwd("/Users/andrespasinetti/Bayesian-Project-4")

source("dags/propose_DAG.R")
source("dags/operation.R")
source("dags/acceptreject_DAG.R")
source("dags/new_bcdag.R")
source("dags/update_DAG.R")

### Plot 3d
#install.packages(c("rgl", "car"))

# PARAMETERS
C <- 3 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
q <- 3
N_SAMPLES <- 700 # number of samples
niter <- 200 # number of Gibbs Sampling iterations

w_real <- gtools::rdirichlet(1, alpha_0)
DAG_real <- array(dim = c(q, q, C))
Omega_real <- array(dim = c(q, q, C))
for (c in 1:C) {    
    DAG_real[, , c] <- rDAG(q = q, w = 0.5)
    L <- matrix(runif(n = q*(q), min = -10, max = 10), q, q)     ### va bene mettere questi min e max? 
    L <- L * DAG_real[, , c] 
    diag(L) <- 1
    D <- diag(1, q)
    Omega_real[, , c] <- L%*%solve(D)%*%t(L)       
}


rmix <- function(w_real, Omega_real) {
    z_real <- sample(1:C, prob = w_real, size = N_SAMPLES, replace = TRUE)
    x <- matrix(, nrow = N_SAMPLES, ncol = q)
    for (i in 1:N_SAMPLES) {
        x[i, ] <- mvtnorm::rmvnorm(1, sigma=solve(Omega_real[, ,z_real[i]]))
    }
  
  return(list(x, z_real))
}

mix <- rmix(w_real, Omega_real)
x <- mix[1]
z_real <- mix[2]

x <- data.frame(x)
z_real <- unlist(z_real)
if (q==2) {
  x11()
  plot(x, type="p", col=z_real, pch=20)
} else if (q==3) {
  open3d()
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= z_real, surface = FALSE)

}



# Full conditionals
sample_w <- function(alpha_0, N) {
  alpha_new <- alpha_0 + N
  w <- gtools::rdirichlet(1, alpha_new)
  return(w)
}

sample_z <- function(mu, Omega, w, x) {
  z <- array(dim = N_SAMPLES)
  for (i in 1:N_SAMPLES) {
    prob <- c()
    summation <- 0
    for (c in 1:C) {
      summation <- summation + w[c] * dmvnorm(x[i, ], mu[c, ], solve(Omega[, , c]))
    }
    # avoid division by 0
    if (summation != 0) {
      for (c in 1:C) {
        prob <- c(prob, w[c] * dmvnorm(x[i, ], mu[c, ], solve(Omega[, , c])) / summation)
      }
    } else {
      prob <- runif(n = C, min = 0, max = 1)
    }

    z[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(z)
}


# Gibbs Sampler
gibbs <- function(x, niter, C, alpha_0) {
  x <- data.matrix(x)
  mu <- matrix(0, nrow = C, ncol = q)
  # Hyperparameters
  w_dags <- 0.5 # probability of adjacency in dag

  # Initialize DAG
  DAG <- array(dim = c(q, q, C))
  Omega <- array(dim = c(q, q, C))
  for (c in 1:C) {    
    DAG[, , c] <- rDAG(q = q, w = w_dags)
    L <- matrix(runif(n = q*(q), min = -10, max = 10), q, q)     ### va bene mettere questi min e max? 
    L <- L * DAG[, , c] 
    diag(L) <- 1
    D <- diag(1, q)
    Omega[, , c] <- L%*%solve(D)%*%t(L)       
  }
  
 
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) * 0.01
 
  # Initialize cluster allocation
  z <- array(dim = N_SAMPLES)
  for (i in 1:N_SAMPLES) {
    z[i] <- sample(1:C, prob = w, size = 1)
  }
  
  N <- rep(0, C) # Samples per cluster

  # Save the Markov Chain 
  w_GS <- array(dim = c(C, niter) )
  w_GS[, 1] <- w

  z_GS <- array(dim = c(N_SAMPLES , niter))
  z_GS[, 1] <- z

  DAG_GS <- array(dim = c(q, q, C, niter))
  DAG_GS[, , , 1] <- DAG

  Omega_GS <- array(dim = c(q, q, C, niter))
  Omega_GS[, , , 1] <- Omega

  
  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 5 == 0) {
      cat(i, "/", niter, "\n")
    }

    for (c in 1:C) {
      N[c] <- sum(z == c)
    }

    
    
    w <- sample_w(alpha_0, N)
    z <- sample_z(mu, Omega, w, x)
    out <- update_DAGS(DAG = DAG, data = x, z = z, a = q, U=diag(1, q), w = w_dags)
    DAG <- out$DAG
    Omega <- out$Omega

    w_GS[, i] <- w
    z_GS[, i] <- z
    DAG_GS[, , , i] <- DAG
    Omega_GS[, , , i] <- Omega
  }

  return(list("w_GS" = w_GS, "z_GS" = z_GS, "DAG_GS" = DAG_GS, "Omega_GS" = Omega_GS))
}

mix <- gibbs(x, niter, C, alpha_0)


if (q==2) {
  x11()
  plot(x, type="p", col=mix$z_GS[, niter], pch=20)
} else if (q==3) {
  open3d()
  x <- data.frame(x)
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= mix$z_GS[, niter], surface = FALSE)

}


cat("\nDAG_real:\n")
print(DAG_real)
cat("\nDAG_GS:\n")
print(mix$DAG_GS[, , , niter])
cat("\nOmega_real:\n")
print(Omega_real)
cat("\nOmega:\n")
print(mix$Omega_GS[, , , niter])


#x <- data.frame(x)
