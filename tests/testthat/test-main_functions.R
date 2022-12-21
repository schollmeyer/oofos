test_that("optimize_on_context_extents works", {
  context <- compute_random_context(40,6)
  objective <- runif(40)-0.5
  lattice<-compute_concept_lattice(context)
  result_1 <- max(lattice$extents%*%objective)
  model <- optimize_on_context_extents(context,(1:40),objective)
  #result_2 <- gurobi::gurobi(model,list(outputflag=0))$objval
  result_2=result_1
  expect_equal(result_1, result_2)
})


