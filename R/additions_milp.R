add_two_object_implications <- function(model, X, index1, index2) {
  m <- dim(X)[1]
  n <- dim(X)[2]

  t <- 1
  N <- length(index1)
  mask2 <- (1:model$n_rows) %in% index2
  ans <- array(0, c(choose(N, 2), m + n))
  for (kk in (1:(N - 1))) {
    k <- index1[kk]
    # for(k in index[(1:(N-1))]){


    idx1 <- (X[k, ] == 1)

    for (l in index1[((kk + 1):N)]) {
      idx2 <- (X[l, ] == 1)
      idx <- idx1 & idx2
      idx <<- idx
      if (length(which(idx == 1)) > 0) {
        XX <- matrix(X[, idx], ncol = length(which(idx == 1)))
        i <- which(rowSums(XX) == sum(idx1 & idx2) & mask2)
        i <- setdiff(i, c(k, l))


        if (length(i) > 0) {
          ans[t, i] <- -1 / (length(i))
          ans[t, c(k, l)] <- 1

          t <- t + 1
        }
      }
    }
  }
  t <- t - 1
  model2 <- model
  model2$A <- rbind(model2$A, ans[(1:t), ])
  model2$rhs <- c(model2$rhs, rep(1, t))
  model2$sense <- c(model2$sense, rep("<=", t))

  return(model2)
}


add_two_attribute_implications <- function(model, X) {
  m <- dim(X)[1]
  n <- dim(X)[2]
  t <- 1

  ans <- array(0, c(choose(n, 2), m + n))
  for (k in (1:(n - 1))) {
    idx1 <- (X[, k] == 1)

    for (l in ((k + 1):n)) {
      idx2 <- (X[, l] == 1)
      i <- which(colSums(X[idx1 & idx2, ]) == sum(idx1 & idx2))

      if (length(i) > 2) {
        ans[t, i + m] <- -1 / (length(i) - 2)
        ans[t, c(k + m, l + m)] <- 1

        t <- t + 1
      }
    }
  }
  t <- t - 1
  model2 <- model
  model2$A <- rbind(model2$A, ans[(1:t), ])
  model2$rhs <- c(model2$rhs, rep(1, t))
  model2$sense <- c(model2$sense, rep("<=", t))

  return(model2)
}






optimistic_estimate_of_pairs <- function(X, v) {
  m <- dim(X)[1]
  n <- dim(X)[2]
  i <- (v > 0)
  indexs <- (1:m)
  ans <- array(0, c(n, n))
  for (k in (1:n)) {
    idx1 <- (X[, k] == 1)

    for (l in (k:n)) {
      idx2 <- (X[, l] == 1)
      ans[k, l] <- sum(i * v * (idx1 == 1 & idx2 == 1))
      if (k != l) {
        ans[l, k] <- ans[k, l]
      }
    }
  }

  return(ans)
}



add_sos_constraints <- function(model, v_max) {
  m <- dim(model$context)[1]
  n <- dim(model$context)[2]
  O <- optimistic_estimate_of_pairs(model$context, model$obj[(1:m)])
  I <- (O <= v_max)
  idx <- (which(diag(I)))
  diag(I) <- 0
  graph <- igraph::graph_from_adjacency_matrix((1-I))##igraph::as((1 - I), "graphNEL")
  graph <- igraph::as_graphnel(graph)
  coloring <- RBGL::sequential.vertex.coloring(graph)
  K <- coloring[[1]]
  t <- 1
  sos <- list()
  for (k in (1:K)) {
    indexs <- which(coloring[[2]] == k - 1)
    if (length(indexs) > 1) {
      sos[[t]] <- list(index = indexs + m, type = 1, weight = rep(1, length(indexs)))
      t <- t + 1
    }
  }
  model2 <- model
  model2$sos <- sos
  if (length(idx) >= 1) {
    model2$ub[idx + m] <- 0
  }
  return(model2)
}


add_lazy_constraints <- function(model, v_max, lazy = 1, eps = 0) {
  m <- dim(model$context)[1]
  n <- dim(model$context)[2]
  O <- optimistic_estimate_of_pairs(model$context, model$obj[(1:m)])
  I <- (O <= v_max)

  graph <- igraph::as(1 - I, "graphNEL")
  coloring <- RBGL::sequential.vertex.coloring(graph)
  K <- coloring[[1]]
  t <- 1
  A <- array(as.logical(0), c(K, m + n))
  for (k in (1:K)) {
    indexs <- which(coloring[[2]] == k - 1)
    if (length(indexs) > 1) {
      A[t, indexs + m] <- 1


      t <- t + 1
    }
  }
  t <- t - 1
  model2 <- model
  model2$A <- rbind(model2$A, A[(1:t), ])
  model2$rhs <- c(model$rhs, rep(1 + eps, t))
  model2$sense <- c(model2$sense, rep("<=", t))
  model2$lazy <- c(rep(0, dim(model$A)[1]), rep(lazy, t))


  return(model2)
}




optimistic_estimate_of_pairs2 <- function(X, v, timelimit) {
  M <- extent.opt(X, which(v > 0), v)

  m <- dim(X)[1]
  n <- dim(X)[2]
  i <- (v > 0)
  indexs <- (1:m)
  ans <- array(0, c(n, n))
  for (k in (1:n)) {
    print(k)

    idx1 <- (X[, k] == 1)

    for (l in (k:n)) {
      MM <- M
      MM$lb[k + m] <- 1
      MM$lb[l + m] <- 1


      ans[k, l] <- gurobi::gurobi(MM, list(presolve = 0, outputflag = 0, timelimit = timelimit)$objbound)
      if (k != l) {
        ans[l, k] <- ans[k, l]
      }
    }
  }

  return(ans)
}





extent_opt_feasible <- function(X, gen.index, v, binary.variables = "afap", C, add.sos.constraints) {
  temp <- extent.opt(X = X, gen.index = gen.index, v = v, binary.variables = binary.variables)
  if (add_sos_constraints) {
    temp <- add_sos_constraints(temp, C)
  }
  temp$A <- rbind(temp$A, matrix(temp$obj, nrow = 1))
  temp$sense <- c(temp$sense, ">=")
  temp$rhs <- c(temp$rhs, C)
  return(temp)
}
