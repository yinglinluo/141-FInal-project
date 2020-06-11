

test_that("blblm works", {
  fitlm <- blblm(mpg ~ wt * disp * drat, data= mtcars, m = 3, B = 100)
  expect_s3_class(fitlm, "blblm")
  colm <- coef(fitlm)
  expect_equal(length(colm), 8)
})

