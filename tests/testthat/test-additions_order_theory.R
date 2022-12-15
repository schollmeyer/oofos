test_that("compute_width works", {

  example_posets <- compute_example_posets(8)
  expect_equal(compute_width(example_posets$chain)$width,1)
  expect_equal(compute_width(example_posets$maximum_edge_poset)$width,4)
  expect_equal(compute_width(example_posets$antichain)$width,8)
  expect_equal(compute_width(example_posets$two_dimensional_grid)$width,8)
  expect_equal(compute_width(example_posets$powerset_order)$width,70)

})


test_that("compute_pseudoreduction works", {

  for( n in c(5,10,15,20)){
    prob <- runif(1)
    incidence <- compute_random_context(n,n,prob=prob)
    result_1 <- compute_transitive_hull(compute_pseudoreduction(incidence))
    result_2 <- compute_transitive_hull(incidence)
    expect_equal(all(result_1==result_2),TRUE)
  }
})
