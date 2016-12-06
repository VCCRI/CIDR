context("cpp_dist_weighted - dissimilarity matrix using weighted imputation")

test_that("distance calculation works", {
  N <- 2
  a <- 1.5
  b <- 5.0
  
  # Both dropout candidates
  dist <- array(0, dim=c(N, N))
  dropouts <- array(1, dim=c(1, N))
  data <- array(c(1.2, 0.9), dim=c(1,N))
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(0, 0)), dim=c(2,2)))
  
  # Neither dropouts, same counts
  dist <- array(0, dim=c(N, N))
  dropouts <- array(0, dim=c(1, N))
  data <- array(c(3.5, 3.5), dim=c(1,N))
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(0, 0)), dim=c(2,2)))
  
  # Neither dropouts, different counts
  dist <- array(0, dim=c(N, N))
  dropouts <- array(0, dim=c(1, N))
  data <- array(c(3.5, 3.6), dim=c(1,N))
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(0.01, 0)), dim=c(2,2)))
  
  # #1 dropout, #2 not dropout (imputation occurs)
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(1, 0), dim=c(1, N))
  data <- array(c(0.5, 5.0), dim=c(1,N))
  # P(5.0) = 0.5, imputed value = 2.75
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(5.0625, 0)), dim=c(2,2)))
  
  # #1 not dropout, #2 dropout (imputation occurs)
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(0, 1), dim=c(1, N))
  data <- array(c(3.5, 0.5), dim=c(1,N))
  # P(3.5) = 0.9046505351, imputed value = 3.2139516053
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(0.081823684110447, 0)), dim=c(2,2)))
  
  # 3 features
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(c(0, 0, 0),
                      c(1, 1, 0)), dim=c(3, N))
  data <- array(c(c(5.0, 5.0, 4.5),
                  c(0.5, 1.0, 2.5)), dim=c(3,N))
  cpp_dist_weighted(dist, dropouts, data, N, a, b)
  expect_equal(dist, array(c(c(0, 0), c(13.0625, 0)), dim=c(2,2)))
})
