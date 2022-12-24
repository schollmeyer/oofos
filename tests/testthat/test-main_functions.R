test_that("optimize_on_context_extents works", {
  context <- compute_random_context(40,6)
  objective <- sample(c(0,1),size=nrow(context),replace=TRUE)
  objective <- compute_objective(data.frame(objective), target="objective",
                                 target_class=0)
  lattice<-compute_concept_lattice(context)
  result_1 <- max(lattice$extents%*%objective)
  model <- optimize_on_context_extents(context,(1:40),objective)
  result_2 <- gurobi::gurobi(model,list(outputflag=0))
  #result_2=result_1
  quality <- compute_quality(model,result_2)
  expect_equal(result_1, result_2$objval)
})


