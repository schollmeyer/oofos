test_that("compute_example_contexts works", {
  example_contexts <- compute_example_contexts()
  expect_equal(length(example_contexts),3)
})
