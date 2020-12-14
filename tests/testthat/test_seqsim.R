library(gwpcR)
library(data.table)
context('Sequencing simulation')

test_that("seqsim", {
  expect_length(seqsim(c(1), reads.target=10, efficiency=0), 1)
  expect_length(seqsim(c(1,1), reads.target=10, efficiency=0), 2)
  expect_length(seqsim(c(1,2,1), reads.target=10, efficiency=0), 3)
  expect_length(seqsim(c(1,1,2,1,2,3), reads.target=10, efficiency=0), 6)

  # Check that the output order is not shuffled
  r <- seqsim(c(1,1,1e7,2,1,2,3), reads.target=10, efficiency=0)
  expect_length(r, 7)
  expect_equal(r[3], sum(r))
})
