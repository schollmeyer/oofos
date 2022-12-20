test_that("compute_extent_vc_dimension works", {

  context <- diag(rep(1,10))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,2)

  context <- 1-diag(rep(1,10))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,10)

  context <- upper.tri(diag(1,20))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,1)

  context <- cbind(upper.tri(diag(1,3)),lower.tri(diag(1,3)))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,2)



})
