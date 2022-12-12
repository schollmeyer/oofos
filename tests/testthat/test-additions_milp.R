test_that("add_two_object_implications  works", {
  n_obj <- 30
  context <- compute_random_context(n_obj,7,p=0.8,seed=NULL)
  objective <- runif(n_obj)-0.5
  model <- optimize_on_context_extents(context,(1:n_obj),objective)
  result_1 <- gurobi(model,list(outputflag=0))$objval
  model_2 <- add_two_object_implications(model,context,(1:n_obj),(1:n_obj))
  result_2 <- gurobi(model_2,list(outputflag=0))$objval
  expect_equal(result_1,result_2)
})
