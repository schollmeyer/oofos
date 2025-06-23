
compute_maximal_elements <- function(relation_mat) {
  mat <- relation_mat
  diag(mat) <- 1
  mask <- !duplicated(mat)&!duplicated(t(mat))
  index <- which(mask)
  mat <- mat[index,index]
  diag(mat) <- 0
  return(index[which(rowSums(mat)==0)])
}

compute_cover_relation <- function(incidence) {
  # computes the cover relation of a poset
  # If the relation incidence is transitive and antysmmetric bot not reflexive,
  # then the cover relation of the reflexive hull of incidence is computed.
  diag(incidence) <- 1
  n_rows <- nrow(incidence)
  incidence <- as.logical(incidence - compute_relation_product(
    incidence - diag(rep(1, n_rows)), incidence - diag(rep(1, n_rows))))
  dim(incidence) <- c(n_rows, n_rows)
  diag(incidence) <- 0
  return(incidence)
}

compute_cover_ordered_quoset <- function(incidence) {
  # helper function that computes a cover relation of a quasiordered set
  # where the incidence matrix has to be already ordered in increasing
  # columnsums. This function is used by the function 'compute_pseudoreduction'
  incidence_1 <- incidence * lower.tri(incidence)
  incidence_2 <- incidence * upper.tri(incidence)
  result <- pmax(compute_cover_relation(incidence_1),
                 compute_cover_relation(incidence_2))
  return(result)
}

# TODO : checken ob das nicht doch eine MINIMALE /ECHTE transitive reduktion
# ist
compute_pseudoreduction <- function(incidence) {
  # computes a transitive pseudoreduction of an arbitrary homogeneous
  # relation (result is not necessarly minimal w.r.t. set inclusion and
  # therefore is not a proper transitive reduction, cf., also
  # https://en.wikipedia.org/wiki/Transitive_reduction)
  diagonal_elements <- diag(incidence)
  incidence <- compute_transitive_hull(incidence)
  o <- order(colSums(incidence))
  result <- compute_cover_ordered_quoset(incidence[o, o])
  p <- Matrix::invPerm(o)
  result <- result[p, p]
  diag(result) <- diagonal_elements
  return(result)
}


compute_example_posets <- function(n,powerset_order=TRUE) {
  antichain <- diag(rep(1, n))
  chain <- upper.tri(antichain)
  diag(chain) <- 1
  maximum_edge_poset <- array(0, c(n, n))
  m <- ceiling(n / 2)
  maximum_edge_poset[(1:m), -(1:m)] <- 1
  diag(maximum_edge_poset) <- 1

  q_n <- diag(rep(1, n))
  r_n <- 1 - q_n
  s_n <- lower.tri(q_n)
  diag(s_n) <- 1
  t_n <- upper.tri(q_n)
  diag(t_n) <- 1
  p_n <- rbind(q_n, r_n, s_n[-c(1, n - 1, n), ], t_n[-c(1, 2, n), ])
  kellys_poset <- compute_incidence(p_n)
  # TODO : checken ob richtig
  # literature: David Kelly: On the dimension of partially ordered sets.
  # Discrete Mathematics 35 (1981) 135-156

  two_dimensional_grid <- compute_incidence(gtools::permutations(n, 2,
    repeats.allowed = TRUE
  ))
  if(powerset_order==TRUE){
  powerset_order <- compute_incidence(gtools::permutations(2, n,
    repeats.allowed = TRUE
  ) - 1)}
  else{powerset_order <- NULL}
  # interval_order
  upper <- (1:n)
  lower <- upper - stats::runif(n) * 4.01
  interval_order <- array(0, c(n, n))
  for (k in (1:n)) {
    for (l in (1:n)) {
      interval_order[k, l] <- (upper[k] <= lower[l])
    }
  }
  p5 <- ddandrda::compute_all_partial_orders(5,list=FALSE,complemented=FALSE)
  index <- sample(seq_len(nrow(p5)),size=1)
  random_p_order_1 <- p5[index,]; dim(random_p_order_1) <- c(5,5)
  index <- sample(seq_len(nrow(p5)),size=1)
  random_p_order_2 <- p5[index,]; dim(random_p_order_2) <- c(5,5)
  index <- sample(seq_len(nrow(p5)),size=1)
  random_p_order_3 <- p5[index,]; dim(random_p_order_3) <- c(5,5)


  return(list(
    chain = chain, antichain = antichain,
    maximum_edge_poset = maximum_edge_poset, kellys_poset = kellys_poset,
    two_dimensional_grid = two_dimensional_grid,
    powerset_order = powerset_order, interval_order = interval_order,
    random_p_order_1 = random_p_order_1, random_p_order_2 = random_p_order_2,
    random_p_order_3 = random_p_order_1
  ))
}


compute_incidence <- function(X) {
  # erzeugt Inzidenzmatrix einer gegebenen Datentabelle (Zeilen entsprechen
  # statistischen Einheiten und Spalten entsprechen Auspraegungen verschiedener
  # Dimensionen. statistische Einheit x ist kleinergleich statistische Einheit y
  # iff x_i <= y_i fuer jede Dimension i)
  m <- dim(X)[1]
  ans <- matrix(FALSE, ncol = m, nrow = m)
  for (k in (1:m)) {
    for (l in (1:m)) {
      ans[k, l] <- all(X[k, ] <= X[l, ])
    }
  }
  return(ans)
}

# TODO MILP version von width plus beides geneinander testen. Achtung: bei
# nicht-posets kann ergebnis unterschiedlich sein.

# TODO :? principalfilter VC dim (ist bei distrib verb identisch zu odim,...)

#' Compute the width of a partiallyordered set
#'
#' @description 'compute_width' computes the width of partially ordered set,
#' i.e., the maximal cardinality of an antichain. (An antichain is a set of
#' pairwise incomparable elements.)
#' For computing the width is done with a maximum matching formulation, see
#' \url{https://mathoverflow.net/questions/189161/fastest-algorithm-to-compute-the-width-of-a-poset}
#' If 'incidence' is an arbitrary homogeneous relation, then first the reflexive
#' hull is computed and then the width of the quotient order ( i.e., incidence
#' factorized over incidence cap incidence^T )
#' @param incidence  square incidence matrix with 0-1 entries
#' @param quotient_order logical If true, before computing the width,
#' the quotient relation of the relation incidence is computet
#' @return a list with entry width that gives the width of the poset.
#'
#' @export
compute_width <- function(incidence,quotient_order=TRUE) {
  # TODO : Antiketten auch berechnen und ausgeben (Achtung mit Faktorisierung)
  if (is.null(incidence)) {
    return(list(width = 0))
  }
  n_rows <- nrow(incidence)
  if (n_rows != ncol(incidence)) {
    print("No square matrix")
    return(NULL)
  }
  if (any(!incidence %in% c(0, 1))) {
    print("Not all entries of the incidence matrix are in{0,1}")
    return(NULL)
  }
  incidence <- compute_transitive_hull(incidence)
  diag(incidence) <- 1
  if(quotient_order){incidence <- compute_quotient_order(incidence)}
  n_rows <- nrow(incidence)
  if (n_rows == 0) {
    return(list(width = 0))
  }
  if (n_rows == 1) {
    return(list(width = 1))
  }
  incidence <<- incidence

  diag(incidence) <- 0
  graph_incidence <- rbind(
    cbind(0 * incidence, incidence),
    cbind(0 * incidence, 0 * incidence)
  )
  colnames(graph_incidence) <- as.character(seq_len(2*n_rows))
  rownames(graph_incidence) <- as.character(seq_len(2*n_rows))
  graph <- igraph::graph_from_adjacency_matrix(graph_incidence)
  igraph::V(graph)$type <- c(rep(FALSE, n_rows), rep(TRUE, n_rows))
  ans <- igraph::max_bipartite_match(graph)
  antichain_indexs <- seq_len(n_rows)

  #na_indexs <- which(is.na(ans$matching[seq_len(n_rows)])) #- n_rows
  #antichain_indexs <- antichain_indexs[-na_indexs]
  #antichain_indexs <- setdiff(antichain_indexs,ans$matching[seq_len(n_rows)])
  #antichain_indexs <- antichain_indexs[which(antichain_indexs %in% seq_len(n_rows))]
  return(list(width = n_rows - ans$matching_size,antichain_indexs=antichain_indexs,matching_result=ans,graph_incidence=graph_incidence))
}


compute_width_milp <- function(incidence){
  n_rows <- nrow(incidence)
  diag(incidence) <- 0
  model <- list(A=array(0,c(n_rows,n_rows)), rhs=rep(1,n_rows),
                sense=rep("<=",n_rows),vtypes=rep("B",n_rows),
                obj=rep(1,n_rows),modelsense="max")
  for(k in seq_len(n_rows)){
    model$A[k,which(incidence[k,]==1)]=1
    model$A[k,k] <- sum(incidence[k,])
    model$rhs[k] <- sum(incidence[k,])
  }
  result <- gurobi::gurobi(model,list(outputflag=0))
return(list(width=result$objval,antichain=result$x))}




### # This function (with slight modifications) is taken from the R package
# ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

## TODO : harmonize with latest ddandrda version


compute_relation_product <- function(x, y) {
  # This function is based on code that is part of the master thesis
  # Blocher, Hannah (2020). A Supremum Type Test Statistic for  Distributional
  # Tests: An Analysis in the Context of VC Theory, Partially Spatial Data and
  # Regularization. LMU Munich

  # @X (matrix): Represents a graph with edges (weight one) and knots
  # @Y (matrix): Represents a graph with edges (weight one)and knots
  # Return (matrix): Represents a graph after two steps which are defined by
  #                 X and Y

  # Input check
  if (dim(y)[1] != dim(x)[1]) {
    print("dimension missmatch!")
    stop
  }

  number_row_x <- dim(x)[1]
  number_col_x <- dim(x)[2]
  number_col_y <- dim(y)[2]

  product_result <- array(0, c(number_row_x, number_col_y))

  for (k in (1:number_col_x)) {
    # X[.,k] edge between . to k
    # Y[k,.] edge from k to .

    # computes if there is a path which uses knot k as intermediate step
    product_result[which(x[, k] == 1), which(y[k, ] == 1)] <- 1
  }

  return(product_result)
}



#' Compute the transitive hull of a partial order
#'
#' @description
#'
#' This function is based on code that is part of the master thesis
#' Blocher, Hannah (2020). A Supremum Type Test Statistic for  Distributional
#' Tests: An Analysis in the Context of VC Theory, Partially Spatial Data and
#' Regularization. LMU Munich
#'
#'
#' compute_transitive_hull' returns a 0-1-matrix which represents the
#' order-pairs given by the transitivity property of a partial order
#'
#' @param relation_mat Rrepresents a relation matrix. Note that
#' has to be a squared matrix.
#'
#' @return The transitive hull of the relation matrix relation_mat
#'
#' @examples
#' relation_mat_input <- matrix(0, nrow = 5, ncol = 5)
#' relation_mat_input[1, 3] <- 1
#' relation_mat_input[2, 1] <- 1
#' relation_mat_input[4, 3] <- 1
#' oofos:::compute_transitive_hull(relation_mat_input)
#'
compute_transitive_hull <- function(relation_mat) {
  # @relation_mat (sqared matrix): represents a relation matrix
  # Return (squared matrix): the transitive hull of the relation matrix
  # relation_mat

  number_obj <- dim(relation_mat)[1]

  old_matrix <- array(0, c(number_obj, number_obj))
  next_matrix <- relation_mat
  diag(next_matrix) <- 1 ## TODO  In ddandrda auch korrigieren
  diag(relation_mat) <- 1 ## TODO ...
  transitive_hull <- relation_mat

  # each while-loop computes the next step of the path given by relation_mat
  while (any(old_matrix != next_matrix)) {
    old_matrix <- next_matrix
    # this computes the next step. In other words in the first loop it computes
    # all edges which can be obtained by the combination of twp edges in
    # relation_mat
    next_matrix <- compute_relation_product(old_matrix, relation_mat)

    # contains all paths which can be done in maximal number of loop iteration
    # of steps
    transitive_hull <- transitive_hull + next_matrix
  }

  # more than one possible path --> than the number of paths was computed by the
  # while-loop
  # --> set to 1
  index_non_zero <- which(transitive_hull > 0)
  transitive_hull[index_non_zero] <- 1
  diag(transitive_hull) <- diag(relation_mat)

  return(transitive_hull)
}



#' Compute the order dimension of a partially ordered set
#'
#' @description 'compute_order-dimension' computes the order dimension of a partially ordered set,
#' i.e., the minimal k such that the poset can be represented as the intersection of k linear orderstichain.
#'
#' @param poset square incidence matrix with 0-1 entries that describes the poset
#' @param poset_context For the computation of the order dimension, the context with all partial orders
#' as intents is needed. More precisely, one needs a formal context where every object is a linear order L and
#' every attribute is a pair (a,b) and L I (a,b) iff (a,b) in L. Note that the empty intersection of objects
#' gives the all relation, and the all relation is not a partial order.
#' If poset_context is not supplied, then it will be computed beforehand
#' @return the order dimension of the poset.
#'
#' @export
compute_order_dimension <- function(poset, poset_context=NULL){
  # the context with all partial orders as intents is needed. If it is not supplied, then it will be computed beforehand
  if(is.null(poset_context)){
    n_items <- nrow(poset)
    perms <- gtools::permutations(n_items,n_items)
    poset_context <- ddandrda:::compute_ranking_scaling(perms, remove_full_columns = FALSE, complemented = FALSE)
  }
  # now we compute the formal concept (E,I) with extent E as the set of all linear extensions of the considered poset 'poset' and I as the set of all pairs in 'poset'
  E <- rep(0,nrow(poset_context))
  for(k in seq_len(nrow(poset_context))){
    # if the linear order 'context[k,]' extents is a linear extension of 'poset' then it is in the extent E (i.e., E[k]=1)
    if(all(poset_context[k,] >= as.vector(poset))){E[k] <-1}
  }
  # additionally compute the intent I by applying the derivation operator psi to E
  I <- ddandrda:::calculate_psi(E,poset_context)
  # computes a gurobi modelmodel that itself computes a minimal k such that there exists an object set F of size k that generates E (i.e. gamma(F)=E, where gamma is the closure operator phi \circ psi w.r.t. the context 'context')
  gurobi_model <- min_k_obj_generated(E,I,poset_context)
  #Now the actual gurobi model is computed
  result <- gurobi::gurobi(gurobi_model)
  # the minimal value of the binarry linear program is the minimal k such that the given partial order 'poset' cann be represented as the intersection of k linear orders (concretely the partial orders context[k,] with E[k]=1)
  return(result$objval)}
