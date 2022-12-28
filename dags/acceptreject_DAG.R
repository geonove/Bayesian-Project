acceptreject_DAG <- function(tXX, n, currentDAG, proposedDAG, node, op.type, a, 
                             U, w, current.opcard, proposed.opcard) 
{
  logprior.ratios <- c(log(w/(1 - w)), log((1 - w)/w), log(1))
  logprior.ratio <- logprior.ratios[op.type]
  logproposal.ratio <- log(current.opcard) - log(proposed.opcard)
  if (op.type != 3) {
    current_lml <- BCDAG::DW_nodelml(node, currentDAG, tXX, n, 
                              a, U)
    proposed_lml <- BCDAG::DW_nodelml(node, proposedDAG, tXX, n, 
                               a, U)
  }
  else {
    current_lml <- BCDAG::DW_nodelml(node[1], currentDAG, tXX, 
                              n, a, U) + DW_nodelml(node[2], currentDAG, tXX, 
                                                    n, a, U)
    proposed_lml <- BCDAG::DW_nodelml(node[1], proposedDAG, tXX, 
                               n, a, U) + DW_nodelml(node[2], proposedDAG, tXX, 
                                                     n, a, U)
  }
  acp.ratio <- min(0, proposed_lml - current_lml + logprior.ratio + 
                     logproposal.ratio)
  is.accepted <- log(stats::runif(1)) < acp.ratio
  return(is.accepted)
}

