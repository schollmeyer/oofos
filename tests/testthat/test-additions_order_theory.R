test_that("compute_width works", {

  example_posets <- compute_example_posets(8)
  expect_equal(compute_width(example_posets$chain)$width,1)
  expect_equal(compute_width(example_posets$maximum_edge_poset)$width,4)
  expect_equal(compute_width(example_posets$antichain)$width,8)
  expect_equal(compute_width(example_posets$two_dimensional_grid)$width,8)
  expect_equal(compute_width(example_posets$powerset_order)$width,70)

})
