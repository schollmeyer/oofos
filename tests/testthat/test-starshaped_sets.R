test_that("properties betweenness works", {
  context <- compute_random_context(10,10)
  context <- cbind(context,1-context)
  context <- unique(context)
  betweenness <- get_non_stylized_betweenness(context)
  expect_equal(check_cond_antisymmetry(betweenness),TRUE)
  expect_equal(check_conditional_reflexivity(betweenness),TRUE)
  expect_equal(check_outer_symmetry(betweenness),TRUE)

})

test_that("get_stylized_btweenness works", {
  context <- compute_random_context(10,8)
  betweenness_1 <- get_whole_stylized_betweenness(context)
  betweenness_2 <- get_stylized_betweenness(context)
  betweenness_3 <- get_non_stylized_betweenness(context)
  expect_equal(all(betweenness_1==betweenness_2),TRUE)
  expect_equal(all(betweenness_3 == 1*(betweenness_1 >= 1)),TRUE)
})

test_that("fit_ks_distribution works", {
x  <- rep(0,1000)

for(k in seq_len(1000)){

  x[k] <- (ks.test(rnorm(100),pnorm)$statistic)
}
suppressWarnings( result <- fit_ks_distribution(x,FALSE))

expect_equal(result$value<=0.0001
,TRUE)

})

test_that("compute_starshaped_distr_test works", {
  betweenness <- get_betweenness_from_poset(compute_example_posets(15,powerset_order=FALSE)$two_dimensional_grid)
  n_rows <- nrow(betweenness)
  widths <-compute_widths(betweenness)


  index <- rep(0,n_rows)
  center_id <- sample(seq_len(n_rows),size=1)
  for(k in (1:8)){
     boundary_id <- sample(seq_len(n_rows),size =1)
     index[which(betweenness[center_id,,boundary_id]==1)] <-1
  }
  probs <- rep(1,n_rows)
  probs[which(index==1)] <- 10
  probs <- probs/sum(probs)

  positives <- sample(seq_len(n_rows),size=n_rows,replace=TRUE,prob=probs)
  response <- rep(0,n_rows)
  response[positives] <- 1

  objective <- compute_objective(data.frame(response), target="response",
                                 target_class=1)


  discovery <- discover_starshaped_subgroups(betweenness,objective,complexity_measure=compute_width,complexity_control=Inf)
  test <- compute_starshaped_distr_test(discovery,n_rep=100)

  expect_equal((abs(test$p_value-test$p_value_parametric)<=0.05),TRUE)



})




test_that("discover_starshaped_subgroups works", {
  betweenness <- get_betweenness_from_poset(
    compute_example_posets(6)$two_dimensional_grid
  )
  objective <- stats::rnorm(nrow(betweenness))
  result <- discover_starshaped_subgroups(betweenness, objective,complexity_measure=compute_width, complexity_control= Inf)
  expect_equal(check_if_starshaped(result$star, betweenness), TRUE)
  result_2 <- discover_starshaped_subgroups(betweenness, objective,complexity_measure=compute_width,complexity_control= 8)
  result$objval > result_2$objval


  ternary_relation <- runif(20^3)
  dim(ternary_relation) <- rep(20, 3)
  objective <- rnorm(20)
  result_1 <- discover_starshaped_subgroups(ternary_relation, objective, complexity_measure=compute_width,complexity_control=Inf)
  result_2 <- discover_starshaped_subgroups(ternary_relation, objective, complexity_measure=compute_width,complexity_control=8)
  expect_equal(result_1$objval >= result_2$objval, TRUE)
  expect_equal(result_1$objval > result_2$objval, TRUE)


  n <- 3
  betweenness <- get_betweenness_from_poset(
    compute_example_posets(n)$two_dimensional_grid
  )
  objective <- stats::rnorm(nrow(betweenness))
  result <- discover_starshaped_subgroups(betweenness, objective, complexity_measure=compute_width,complexity_control=Inf)
  sets <- gtools::permutations(2, n^2, repeats.allowed = TRUE) - 1
  objval <- -Inf
  for (k in seq_len(nrow(sets))) {
    if (check_if_starshaped(sets[k, ], betweenness)) {
      objval <- max(objval, sum(objective * sets[k, ]))
    }
  }

  expect_equal(result$objval, objval)
  expect_equal(
    discover_starshaped_subgroups_recompute(result, result$obj),
    result$objval
  )
  result_h0 <- discover_starshaped_subgroups_h0(result)


  ternary_relation <- get_random_ternary_relation(12, prob = 0.5)
  objective <- rnorm(12)
  result <- discover_starshaped_subgroups(ternary_relation, objective, complexity_measure=compute_width,complexity_control=Inf)
  sets <- gtools::permutations(2, 12, repeats.allowed = TRUE) - 1
  objval <- -Inf
  for (k in seq_len(nrow(sets))) {
    if (check_if_starshaped(sets[k, ], ternary_relation)) {
      objval <- max(objval, sum(objective * sets[k, ]))
    }
  }

  result$objval == objval
  expect_equal(result$objval, objval)
})
