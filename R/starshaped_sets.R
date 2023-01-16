##########################################
######## 			     	              ########
####### 						                 #######
######     starshaped sets          ######
#######                            #######
######## 				                  ########
##########################################

check_outer_symmetry <- function(betweenness) {
  n_rows <- nrow(betweenness)
  for (l in seq_len(n_rows)) {
    if (any(betweenness[, l, ] != t(betweenness[, l, ]))) {
      return(FALSE)
    }
  }
  return(TRUE)
}

check_cond_antisymmetry <- function(betweenness) {
  n_rows <- nrow(betweenness)
  for (k in seq_len(n_rows)) {
    mat <- betweenness[k, , ]
    diag(mat) <- 0
    if (any(pmin(mat, t(mat)) == 1)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

check_conditional_reflexivity <- function(betweenness) {
  n_rows <- nrow(betweenness)
  for (k in seq_len(n_rows)) {
    if (any(diag(betweenness[k, , ]) == 0) | any(diag(betweenness[k, , ]) == 0) |
      any(diag(betweenness[, , k]) == 0)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

get_random_ternary_relation <- function(n, prob, reflexive = FALSE,
                                        transitive = FALSE) {
  result <- (stats::runif(n^3) >= prob) * 1
  dim(result) <- rep(n, 3)
  for (k in (1:n)) {
    if (reflexive) {
      diag(result[k, , ]) <- 1
    }
    if (transitive) {
      result[k, , ] <- compute_transitive_hull(result[k, , ])
    }
  }
  return(result)
}

get_betweenness_from_poset <- function(p_order) {
  ## TODO : conditional REFLEXIVITY!!
  n_rows <- nrow(p_order)
  result <- array(0, c(rep(n_rows, 3)))
  for (k in seq_len(n_rows)) {
    for (l in seq_len(n_rows)) {
      for (m in seq_len(n_rows)) {
        result[k, l, m] <- (p_order[k, l] & p_order[l, m]) |
          (p_order[m, l] & p_order[l, k])
      }
    }
  }
  return(result)
}

# compute_examples_starhaped_bet
check_if_starshaped <- function(set, ternary_relation) {
  if (all(set == 0)) {
    return(FALSE)
  }
  n_rows <- nrow(ternary_relation)
  for (k in which(set == 1)) {
    if (check_if_center_point(k, set, ternary_relation)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

check_if_center_point <- function(point, set, ternary_relation) {
  n_rows <- nrow(ternary_relation)
  for (k in which(set == 1)) {
    indexs <- which(ternary_relation[point, , k] == 1)

    if (any(set[indexs] == 0)) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Compute the canoncical (stylized) betweennes relation between objects of a
#' formal context
#'
#' @description 'compute_stylized_betweenness' computes the canonical (stylized)
#' betweennes relation between objects of a formal context: Given three
#' objects g,h,i of a formal context we say that h lies betweeen g and i iff
#'
#'
#' Psi({g}) cap Psi({i}) subseteq Psi({h}).
#'
#'
#' Additionally we say in a stylized
#' manner that h lies betweeen g and i iff
#' Psi({g}) cap Psi({i}) subseteq Psi({h}) is almost true with the exception of
#' some attributes of Psi({g}) cap Psi({i}) that do not belong to Psi({h}) but
#' that are not so 'important' in the sense that there are not too much objects
#' that have all attributes in (Psi({g}) cap Psi({i})) -  Psi({h}) but that have
#' not all attributes in Psi({h}). Concretely the stylized betweennes is
#' quantified by the maximum of the weights 'attribute_weights' associated
#' to these attributes. (The weights could be defined for example as the cloumn
#' means of the underlying formal context).
#'
#' @param g is the object g represented by a 0-1 vector of attributes.
#'
#' @param h is the object h represented by a 0-1 vector of attributes.
#'
#' @param i is the object i represented by a 0-1 vector of attributes.
#' @param context is the underlying context.
#' @param attribute_weights is the weight vector for the attributes.
#' @return a ternary fuzzy relation given by an array of dimension
#' n x n x n where n is the number of rows of the context. Higher values
#' of an entry correspond to a larger degree of betweenness.
#' @export
compute_stylized_betweenness <- function(g, h, i, context, attribute_weights =
                                           colMeans(context)) {
  common_attributes <- which(g == 1 & i == 1)
  if (length(common_attributes) == 0) {
    return(1)
  }
  ans <- 1 - max((1 - h[common_attributes]) *
    attribute_weights[common_attributes])
  return(ans)
}

get_whole_stylized_betweenness <- function(context,
                                           stylized_betweenness_function
                                           = compute_stylized_betweenness,
                                           ...) {
  n_rows <- nrow(context)
  result <- array(0, rep(n_rows, 3))
  pb <- utils::txtProgressBar(min = 0, max = n_rows, initial = 0)
  for (k in (1:n_rows)) {
    utils::setTxtProgressBar(pb, k)
    # if(print_progress){print(c("progress"))}
    for (l in (1:n_rows)) {
      for (m in (1:n_rows)) {
        result[k, l, m] <- stylized_betweenness_function(
          context[k, ], context[l, ],
          context[m, ], context, ...
        )
      }
    }
  }
  close(pb)
  return(result)
}

get_non_stylized_betweenness <- function(context) {
  n_rows <- nrow(context)
  betweenness <- array(0, rep(n_rows, 3))
  pb <- utils::txtProgressBar(min = 0, max = n_rows, initial = 0)
  for (k in seq_len(n_rows)) {
    utils::setTxtProgressBar(pb, k)
    for (m in seq_len(n_rows)) {
      temp <- pmin(context[k, ], context[m, ])
      for (l in seq_len(n_rows)) {
        betweenness[k, l, m] <- all(temp <= context[l, ])
      }
    }
  }

  close(pb)
  return(betweenness)
}

get_stylized_betweenness <- function(context) {
  n_rows <- nrow(context)
  n_cols <- ncol(context)
  col_means <- colMeans(context)

  betweenness <- array(0, rep(n_rows, 3))
  pb <- utils::txtProgressBar(min = 0, max = n_rows, initial = 0)
  for (k in seq_len(n_rows)) {
    utils::setTxtProgressBar(pb, k)
    for (m in seq_len(n_rows)) {
      temp <- rep(0, n_cols)
      i <- which(context[k, ] == 1 & context[m, ] == 1)
      temp[i] <- col_means[i]
      for (l in seq_len(n_rows)) {
        betweenness[k, l, m] <- max(temp - context[l, ])
      }
    }
  }

  close(pb)
  return(1 - betweenness)
}


#' Compute the quotient order of a quasiorder
#'
#' @description 'compute_quotient_order' computes the quotient order of a
#' quasiorder by deleting equivalent 'duplicates' from the given incidence
#' matrix of the quasiorder. if 'incidence' is not a quasiorder, then first the
#' transitive reflexive hull is computet.
#'
#' @param incidence is the incidence matrix of the quasiorder
#'
#' @return the incidence matrix of the quotient order
#'
#' @export
compute_quotient_order <- function(incidence) {
  diag(incidence) <- 1
  incidence <- compute_transitive_hull(incidence)
  index <- which(!duplicated(incidence) & !duplicated(t(incidence)))
  result <- incidence[index, index]
  dim(result) <- rep(length(index), 2)
  return(result)
}





cut_incidence <- function(incidence, width, interval = stats::quantile(
                            unique(as.vector(incidence)), c(0, 1)
                          )) {
  vc <- compute_width(compute_quotient_order(
    compute_transitive_hull(incidence >= interval[2])
  ))$width
  if (vc <= width) {
    return(incidence >= interval[2])
  }
  f <- function(C, incidence) {
    width_2 <- compute_width(
      compute_quotient_order(compute_transitive_hull(
        incidence >= C
      ))
    )$width
    return(width_2 - width)
  }

  ans <- stats::uniroot(f, interval = interval, incidence = incidence)
  return(compute_transitive_hull(incidence >= ans$root))
}


#' Perform a starshaped subgroup discovery
#'
#' @description 'discover_starshaped_subgroups' performs a starshaped subgroup
#' discovery, see Schollmeyer et al. 2023. An objective function (e.g., the
#' Shapiro- Piatetsky quality function) is maximized over the family of all
#' subgroups that are starshaped w.r.t. a 'fuzzy' betweennes relation
#'
#' TODO : explain VC trimming...
#'
#' @param stylized_betweenness is the fuzzy betweennness relation.
#'
#' @param objective is the (linear) objective function to optimize given as an
#' objective vector of the same length as the number of objects.
#'
#' @param local_vc_dimension Is the local VC dimension: Given a centerpoint
#' c, the VC dimension of the subfamily of all starshaped sets with centerpoint
#' c is controlled with this parameter. If local_vc_dimension is set to Inf,
#' then the complexity of the family of all starshaped sets is not reduced at
#' all.
#'
#' @param params is a list with further arguments that are passed to the gurobi
#' optimizer. By default the outputflag of the gurobi solver is set to 0 (this
#' means that not optimization details are printed during the optimization).
#'
#' @return a list with the following entries: TODO
#' @export
discover_starshaped_subgroups <- function(stylized_betweenness, objective,
                                          local_vc_dimension,
                                          params = list(Outputflag = 0)) {
  if (dim(stylized_betweenness)[1] != dim(stylized_betweenness)[2] |
    dim(stylized_betweenness)[1] != dim(stylized_betweenness)[3] |
    dim(stylized_betweenness)[2] != dim(stylized_betweenness)[3]) {
    print("dimension mismatch")
  }
  n_rows <- nrow(stylized_betweenness)
  used_betweenness <- array(0,rep(n_rows,3))
  max_value <- max(stylized_betweenness)
  model <- list(
    modelsense = "max", obj = objective, lb = rep(0, n_rows),
    ub = rep(1, n_rows)
  )
  solutions <- list()
  objvals <- rep(0, n_rows)
  stars <- array(0, c(n_rows, n_rows))
  runtimes <- list()
  models <- list()
  pb <- utils::txtProgressBar(min = 0, max = n_rows, initial = 0)
  for (k in seq_len(n_rows)) { ## quantify over all starcenters
    utils::setTxtProgressBar(pb, k)

    if (local_vc_dimension == Inf) {
      incidence <- (stylized_betweenness[k, , ] >=
        max_value) * 1
    } else {
      incidence <- cut_incidence(
        stylized_betweenness[k, , ],
        local_vc_dimension
      )
    }

    used_betweenness[k,,] <- incidence


    model <- get_model_from_quasiorder(t(incidence))
    if (is.null(model)) {
      model <- list(
        A = matrix(0, nrow = 1, ncol = n_rows), rhs = 1,
        sense = "<="
      )
    }
    model$A <- slam::as.simple_triplet_matrix(model$A)
    model$obj <- objective
    model$lb <- rep(0, n_rows)
    model$ub <- rep(1, n_rows)
    # force centerpoint to be in the set
    model$lb[k] <- 1
    model$modelsense <- "max"

    b <- gurobi::gurobi(model, params = params)
    solutions[[k]] <- b
    runtimes[[k]] <- b$runtime
    objvals[k] <- b$objval
    stars[k, ] <- b$x

    models[[k]] <- model
  }
  i <- which.max(objvals)



  if (local_vc_dimension == Inf) {
    incidence <- (stylized_betweenness[i, , ] >=
      max(stylized_betweenness[k, , ])) * 1
  } else {
    incidence <- cut_incidence(stylized_betweenness[i, , ], local_vc_dimension)
  }


  model <- get_model_from_quasiorder(t(incidence))
  if (is.null(model)) {
    model <- list(A = matrix(0, nrow = 1, ncol = n_rows), rhs = 1, sense = "<=")
  }
  model$obj <- objective
  model$lb <- rep(0, n_rows)
  model$ub <- rep(1, n_rows)
  # force centerpoint to be in the set
  model$lb[k] <- 1
  model$modelsense <- "max"
  b <- gurobi::gurobi(model, params = params)
  close(pb)
  return(list(
    models = models, obj = objective, solutions = solutions,
    objvals = objvals, stars = stars, objval = objvals[i],
    star = stars[i, ], center_id = i,
    fuzzy_incidence = stylized_betweenness[i, , ],
    incidence = incidence, model = model, runtimes=runtimes, used_betweenness =
      used_betweenness
  ))
}


compute_widths <- function(ternary_relation) {
  # TODO : Antiketten auch berechnen und ausgeben (Achtung mit Faktorisierung)
  n_rows <- nrow(ternary_relation)
  result <- rep(0, n_rows)
  for (k in seq_len(n_rows)) {
    result[k] <- compute_width(ternary_relation[k, , ])$width
  }
  return(list(widths = result))
}

# plot_stars <- function(starshaped_result,distance_function){
#
#   i <- which(starshaped_result$star==1)
#
#   j <- oofos::compute_maximal_elements(starshaped_result$incidence[i,i])
#
#
#   i[j]
#   print(i[j])
#
#
#
#
#

# TODO compute_maximal elements?
# }





discover_starshaped_subgroups_recompute <- function(ssd_result, objective) {
  result <- -Inf
  for (k in seq_len(length(ssd_result$models))) {
    model <- ssd_result$models[[k]]
    model$obj <- objective
    result <- max(result, gurobi::gurobi(model, param = list(
      outputflag = 0
    ))$objval)
  }

  return(result)
}


discover_starshaped_subgroups_h0 <- function(ssd_result, params =
                                               list(outputflag = 0)) {
  v <- sample(ssd_result$models[[1]]$obj)

  result <- -Inf
  for (k in seq_len(length(ssd_result$models))) {
    model <- ssd_result$models[[k]]
    model$obj <- v
    result <- max(result, gurobi::gurobi(model, params = params)$objval)
  }

  return(result)
}
#' Compute a statistical significance test for a starshaped subgroup discovery
#'
#' @description 'compute_starshaped_distr_test' performs a statistical
#' significance test to assess if the largest differences of relative
#' frequencies over all starshaped subgroups that was found in a starshaped
#' subgroup discovery is statistically significantly diefferent from zero.
#'
#' @param ssd_result is the result of a starshaped subgroup discovery that was
#' computed with the function 'discover_starshaped_subgroups'.
#' @param n_rep is the number of permutations in the observation-randomization
#' test.
#'
#' @param plot_progress if TRUE (default) then the progress of the
#' observation-randomization test is plot by drawing the empirical distribution
#' function of the currently computed resampled objective values. Additionally,
#' nonparametric and parametric estimates of the p-value based on the currently
#' computed resamples are printed. TODO : beta genauer
#'
#' @return A list with entries 'objvalues': the resampled objective values,
#' 'p_value'  the (nonparametrically) estimated p-value of the
#' observation-randomization test, 'p_value_parametric': a parametric estimate
#' of the p-value based on a beta-approximation of the distribution of the
#' resampled objective values.
#'
#'
#' @export
compute_starshaped_distr_test <- function(ssd_result, n_rep = 1000,
                                          plot_progress = TRUE) {
  if (n_rep <= 2) {
    print("Be serious!")
    return(NULL)
  }
  objvalues <- rep(0, n_rep)
  for (k in seq_len(n_rep)) {
    objvalues[k] <- discover_starshaped_subgroups_h0(ssd_result)
    x <- objvalues[seq_len(k)]
    p_value <- mean(x >= ssd_result$objval)
    if (k > 2) {
      suppressWarnings(fit <- fit_ks_distribution(x))

      p_value_parametric <- stats::pbeta(ssd_result$objval, fit$par[1],
        fit$par[2], fit$par[3],
        lower.tail = FALSE
      )
    }

    if (plot_progress == TRUE & k > 2) {
      plot(stats::ecdf(x),
        do.points = FALSE, col.01line = NULL,
        main = paste(
          "observed value:",
          round(ssd_result$objval, 4),
          "p-palue:", round(p_value, 4), "\n; param. p-value:",
          round(p_value_parametric, 4), "; n:", k, "\nmedian:", round(stats::median(x),4)
        ), verticals = TRUE,
        xlab = "test statistic", ylab = "cdf (black), density (grey)",
        xlim = c(0, 1.05 * max(c(x, ssd_result$objval)))
      )
      graphics::abline(v = ssd_result$objval, col = "darkblue")
      graphics::abline(v = stats::median(x), col = "darkgreen", lty = 2)
      sort_x <- seq(0, 1, length.out = 1000)
      density_parametric <- stats::dbeta(
        sort_x, fit$par[1], fit$par[2],
        fit$par[3]
      )
      graphics::lines(sort_x, density_parametric / max(density_parametric),
        col = "darkgreen"
      )
      cdf_parametric <- stats::pbeta(sort_x, fit$par[1], fit$par[2], fit$par[3])
      graphics::lines(sort_x, cdf_parametric, col = "darkgreen")
      f <- stats::density(x)
      graphics::lines(f$x, f$y / max(f$y), col = "grey")
      # TODO :result_temp <<- list(objvalues = x, p_value = p_value,
      #                           p_value_parametric =
       #   p_value_parametric)
    }
  }
  #rm(result_temp)
  return(list(
    objvalues = objvalues, p_value = p_value, p_value_parametric =
      p_value_parametric
  ))
}

#' Approximate the distribution of a Kolmogorov-Smirnov types test statistic by
#' the (noncentral) beta distribution
#' @description TODO
#' @param objvalues sampled values of a Kolomogorov Sminrov type test statistic
#'
#' @param plot_result If true the empirical distribution function of the
#' 'objvalues'-values together with the cdf of the fitted beta distribution is
#' plotted.
#' description
#'
#'
#' @references Zhang, J., Wu, Y. Beta Approximation to the Distribution of
#' Kolmogorov-Smirnov Statistic. Annals of the Institute of Statistical
#' Mathematics 54, 577–584 (2002). https://doi.org/10.1023/A:1022463111224
#'
fit_ks_distribution <- function(objvalues, plot_result = FALSE) {
  n_objvalues <- length(objvalues)
  cdf <- stats::ecdf(objvalues)
  x <- seq(0, 1, 0.001)
  y <- cdf(x)
  weights <- rep(1, length(x))
  weights[which(y <= 0.5)] <- 0.5
  weights[which(y >= 0.95)] <- 1.5
  weights <- weights / sum(weights)
  loss <- function(par, x, y, weights) {
    y_hat <- stats::pbeta(q = x, shape1 = max(10^-10,par[1]), shape2 = max(10^-10,par[2]), ncp = max(10^-10,par[3]))
    if(any(is.na(y_hat))){print(par);break}
    return(mean(weights * (y - y_hat)^2))
  }
  # TODO : Optimierung besser machen !!!
  result <- stats::optim(
    par = c(1.1, 1.1, 0), fn = loss, x = x, y = y, weights = weights,  method="L-BFGS-B",lower=c(0,0,0),control = list(maxit = 10000000,method="CG",lower=c(0,0,0))
  )
  if (plot_result) {
    y_hat <- stats::pbeta(
      q = x, shape1 = result$par[1], shape2 = result$par[2],
      ncp = result$par[3]
    )
    plot(stats::ecdf(objvalues),
      do.points = FALSE, col.01line = NULL,
      verticals = TRUE
    )
    graphics::lines(x, y_hat, col = "darkblue")
  }
  return(result)
}

get_model_from_quasiorder <- function(quasiorder) {
  # constructs linear program for the optimization over all upsets of
  # a quasiordered set quasiorder

  reduction <- compute_pseudoreduction(quasiorder)
  n_constraints <- sum(reduction)
  n_rows <- dim(reduction)[1]
  constr_mat <- array(as.integer(0), c(n_constraints, n_rows))
  t <- 1
  sense <- rep("<=", n_constraints)
  for (k in (1:(n_rows - 1))) {
    for (l in ((k + 1):n_rows)) {
      if (reduction[k, l] == 1 & reduction[l, k] == 0) {
        constr_mat[t, k] <- 1
        constr_mat[t, l] <- -1
        t <- t + 1
      }

      if (reduction[k, l] == 0 & reduction[l, k] == 1) {
        constr_mat[t, k] <- -1
        constr_mat[t, l] <- 1
        t <- t + 1
      }

      if (reduction[k, l] == 1 & reduction[l, k] == 1) {
        constr_mat[t, k] <- 1
        constr_mat[t, l] <- -1
        sense[t] <- "="
        t <- t + 1
      }
    }
  }
  t <- t - 1

  if (t == 0) {
    return(NULL)
  }
  constr_mat <- constr_mat[(1:t), ]

  if (t == 1) {
    dim(constr_mat) <- c(1, length(constr_mat))
  }

  ans <- list(
    A = constr_mat, rhs = rep(0, t), sense = sense[(1:t)], lb = rep(0, n_rows),
    ub = rep(1, n_rows)
  )
  return(ans)
}




# ###### stilisierte Zwischenrelationen
#
#
# sb1=function(XR,p){
#   m=dim(XR)[1]
#   n=dim(XR)[2]
#   A=array(as.integer(0),c(m,m,m))
#   B=A
#   C=A
#   D=as.matrix(dist(XR,method="minkowski",upper=TRUE,p=p))
#   MAXD=max(D)#Q=quantile(D,QQ)
#
#   for(K in (1:m)){#(1:(m-2))){
#     for(L in (1:m)){#((K+1):(m-1))){
#       for(M in (1:m)){#((L+1):m)){
#         #A[K,L,M]=max(0,D[L,K]-D[M,K],D[L,M]-D[K,M])
#         A[K,L,M]=sum(D[K,M])/(sum(D[K,L]+D[L,M]))
#         B[K,L,M]=(D[K,L]<D[K,M])#*D[K,M]/MAXD#<=Q
#
#         #temp=sign((XR[L,]-XR[K,])*(XR[M,]-XR[K,]))
#         #temp2=abs(XR[L,]-XR[K,]) < abs(XR[M,]-XR[K,])
#         #i=which(temp>=0 &temp2)
#         #A[K,L,M]=length(i)
#         #B[K,L,M]=sum(abs(temp[i]))
#
#       }}}
#   for(K in (1:m)){
#     diag(A[K,,])=1
#     diag(B[K,,])=1}
#   return(list(A=A,B=B))}
#
#
# #######
# #######
# #######
