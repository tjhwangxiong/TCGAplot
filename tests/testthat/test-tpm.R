test_that("get tpm matrix", {
  coad=get_tpm("COAD")
  expect_true(unique(coad$Cancer)=="COAD")
})
