test_that("expand_formula returns a string", {
  expect_true(
    class(expand_formula(~a:b + d, a=c('a1','a2'), b=c('b1','b2'), d='d1'))=="character"
  )
})
