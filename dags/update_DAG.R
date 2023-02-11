update_DAG <- function(DAG, data, a, U, w, fast = FALSE, 
                       collapse = FALSE, verbose = TRUE) {
  
  
  input <- as.list(environment())
  data_check <- sum(is.na(data)) == 0
  a_check <- is.numeric(a) & (length(a) == 1) & (a > ncol(data) - 1)
  w_check <- is.numeric(w) & (length(w) == 1) & (w <= 1) & 
    (w >= 0)
  U_check <- is.numeric(U) & (dim(U)[1] == dim(U)[2]) & (prod(eigen(U)$values) > 
                                                           0) & isSymmetric(U)
  U.data_check <- dim(U)[1] == ncol(data)
  if (data_check == FALSE) {
    stop("Data must not contain NAs")
  }
  if (a_check == FALSE) {
    stop("a must be at least equal to the number of variables")
  }
  if (w_check == FALSE) {
    stop("w must be a number between 0 and 1")
  }
  if (U_check == FALSE) {
    stop("U must be a squared symmetric positive definite matrix")
  }
  if (U.data_check == FALSE) {
    stop("U must be a squared spd matrix with dimensions equal to the number of variables")
  }
  X <- scale(data, scale = FALSE)
  tXX <- t(X) %*% X
  n <- dim(data)[1]
  q <- dim(data)[2]
  DAG_check <- (q == dim(DAG)[1])
  if (DAG_check == FALSE){
    stop("The dimensions of DAG and data do not correspond!")
  }
 # if (save.memory == TRUE) {
  #  Graphs <- vector("double", n.iter)
#    L <- vector("double", n.iter)
 #   D <- vector("double", n.iter)
  #}
  #else {
  #  Graphs <- array(0, dim = c(q, q, n.iter))
  # L <- array(0, dim = c(q, q, n.iter))
  # D <- array(0, dim = c(q, q, n.iter))
  #}
  currentDAG <- DAG
  type = "collapsed"
  prop <- propose_DAG(currentDAG, fast)
  is.accepted <- acceptreject_DAG(tXX, n, currentDAG, 
                                      prop$proposedDAG, prop$op.node, prop$op.type, 
                                      a, U, w, prop$current.opcard, prop$proposed.opcard)
  if (is.accepted == TRUE) {
    currentDAG <- prop$proposedDAG
    }
  Graphs <- currentDAG
  if (collapse == FALSE) {
    type = "complete"
    postparams <- BCDAG::rDAGWishart(1, Graphs, 
                                a + n, U + tXX)
    L <- postparams$L
    D <- postparams$D
    #print(Graphs)
    Sigma <- solve(t(L))%*%(D)%*%solve(L)
  }

if (collapse == FALSE) {
  out <- new_bcdag(list(Graphs = Graphs, L = L, D = D, Sigma = Sigma), 
                   input = input, type = type)
}
else {
  out <- new_bcdag(list(Graphs = Graphs), input = input, 
                   type = type)
}
return(out)
}


### if (!requireNamespace("BiocManager", quietly = TRUE))
###     install.packages("BiocManager")
### BiocManager::install("graph", version = "3.8")
### 
### install.packages("BiocManager")
### BiocManager::install("Rgraphviz")


# Update DAGS
update_DAGS <- function(DAG, data, z, a, U, w) {
  Sigma <- array(dim = c(q, q, C))

  for (c in 1:C) {
    clus_data <- data[z == c, ]
    if (is.null(dim(clus_data))) {
      print("-----")
      clus_data <- array(data[1, ], dim = c(1, q))
    }
    #print(clus_data)


    out <- update_DAG(DAG[, , c], clus_data, a, U, w, fast = TRUE, collapse = FALSE, verbose = TRUE)
    DAG[, , c] <- out$Graphs
    Sigma[, , c] <- out$Sigma
  }
  
  return (list(DAG = DAG, Sigma = Sigma))
}

