test_that("compute_width works", {
  example_posets <- compute_example_posets(8)
  expect_equal(compute_width(example_posets$chain)$width, 1)
  expect_equal(compute_width(example_posets$maximum_edge_poset)$width, 4)
  expect_equal(compute_width(example_posets$antichain)$width, 8)
  expect_equal(compute_width(example_posets$two_dimensional_grid)$width, 8)
  expect_equal(compute_width(example_posets$powerset_order)$width, 70)
})

test_that("compute_width_milp works", {
example_posets <- compute_example_posets(8)
result_1 <- compute_width(example_posets[[6]])
result_2 <- compute_width_milp(example_posets[[6]])
expect_equal(result_1$width,result_2$width)
result_1 <- compute_width(example_posets[[7]])
result_2 <- compute_width_milp(example_posets[[7]])
expect_equal(result_1$width,result_2$width)
result_1 <- compute_width(example_posets[[5]])
result_2 <- compute_width_milp(example_posets[[5]])
expect_equal(result_1$width,result_2$width)

})

test_that("compute_pseudoreduction works", {
  for (n in c(5, 10, 15, 20)) {
    prob <- runif(1)
    incidence <- compute_random_context(n, n, prob = prob)
    result_1 <- compute_transitive_hull(compute_pseudoreduction(incidence))
    result_2 <- compute_transitive_hull(incidence)
    expect_equal(all(result_1 == result_2), TRUE)
  }
})


test_that("compute_widths works", {
  incidence <- compute_random_context(40, 40)
  betweenness <- get_whole_stylized_betweenness(incidence)
  widths1 <- compute_widths(betweenness >= stats::quantile(
    as.vector(betweenness), 0.8
  ))
  widths2 <- compute_widths(betweenness >= stats::quantile(
    as.vector(betweenness), 0.7
  ))
  expect_equal(all(widths1$widths >= widths2$widths), TRUE)
})
