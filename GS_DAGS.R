# PARAMETERS
C <- 4 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
q <- 2
N_SAMPLES <- 200 # number of samples
niter <- 300 # number of Gibbs Sampling iterations

w_real <- gtools::rdirichlet(1, alpha_0)
DAG <- array(dim = c(q, q, C))
Omega <- array(dim = c(q, q, C))
for (i in 1:C) {    
    DAG[, , i] <- rDAG(q = q, w = 0.5)
    L <- matrix(runif(n = q*(q-1), min = 0.1, max = 1), q, q)       ### PERCHE' q*(q-1) e non q*q
    L <- L * DAG[, , i] 
    diag(L) <- 1
    D <- diag(1, q)
    Omega[, , i] <- L%*%solve(D)%*%t(L)       ### PERCHE' SOLVE(D) ???
}


rmix <- function(w_real, Omega) {
    z_real <- sample(1:C, prob = w_real, size = N_SAMPLES, replace = TRUE)
    x <- matrix(, nrow = N_SAMPLES, ncol = q)
    for (i in 1:N_SAMPLES) {
        x[i, ] <- mvtnorm::rmvnorm(1, sigma=Omega[, ,z_real[i]])
    }
  
  return(list(x, z_real))
}

mix <- rmix(w_real, Omega)
x <- mix[1]
z_real <- mix[2]

x <- data.frame(x)
z_real <- unlist(z_real)
x11()
plot(x, type="p", col=z_real, pch=20)