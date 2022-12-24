test_that("starshaped_subgroup_discovery works", {
  betweenness <- get_betweenness_from_p_order(
    compute_example_posets(6)$two_dimensional_grid
  )
  objective <- stats::rnorm(nrow(betweenness))
  result <- starshaped_subgroup_discovery(betweenness, objective, Inf)
  expect_equal(check_if_starshaped(result$star, betweenness), TRUE)
  result_2 <- starshaped_subgroup_discovery(betweenness, objective, 8)
  result$objval > result_2$objval

  ternary_relation <- runif(20^3)
  dim(ternary_relation) <- rep(20, 3)
  objective <- rnorm(20)
  result_1 <- starshaped_subgroup_discovery(ternary_relation, objective, Inf)
  result_2 <- starshaped_subgroup_discovery(ternary_relation, objective, 8)
  expect_equal(result_1$objval > result_2$objval, TRUE)


  n <- 3
  betweenness <- get_betweenness_from_p_order(
    compute_example_posets(n)$two_dimensional_grid
  )
  objective <- stats::rnorm(nrow(betweenness))
  result <- starshaped_subgroup_discovery(betweenness, objective, Inf)
  sets <- gtools::permutations(2, n^2, repeats.allowed = TRUE) - 1
  objval <- -Inf
  for (k in seq_len(nrow(sets))) {
    if (check_if_starshaped(sets[k, ], betweenness)) {
      objval <- max(objval, sum(objective * sets[k, ]))
    }
  }

  expect_equal(result$objval, objval)
  expect_equal(
    starshaped_subgroup_discovery_recompute(result$models, result$obj),
    result$objval
  )
  result_h0 <- starshaped_subgroup_discovery_h0(result$models)


  ternary_relation <- get_random_ternary_relation(12, prob = 0.5)
  objective <- rnorm(12)
  result <- starshaped_subgroup_discovery(ternary_relation, objective, Inf)
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
