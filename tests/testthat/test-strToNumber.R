test_that("strToNumber converts apporpiately", {
  expect_true(is.integer(strToNumber("123L")))
  expect_identical(strToNumber("65L"), 65L)
  expect_true(is.numeric(strToNumber("412")))
  expect_identical(strToNumber("734534"), 734534)
})
