context("cpp_dist - C++ function to create dissimilarity matrix")

test_that("distance calculation works", {
  N <- 2
  
  # Both dropout candidates
  dist <- array(0, dim=c(N, N))
  dropouts <- array(1, dim=c(1, N))
  data <- array(c(1.2, 0.9), dim=c(1,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(0, 0)), dim=c(2,2)))
  
  # Neither dropouts, same counts
  dist <- array(0, dim=c(N, N))
  dropouts <- array(0, dim=c(1, N))
  data <- array(c(3.5, 3.5), dim=c(1,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(0, 0)), dim=c(2,2)))
  
  # Neither dropouts, different counts
  dist <- array(0, dim=c(N, N))
  dropouts <- array(0, dim=c(1, N))
  data <- array(c(3.5, 3.6), dim=c(1,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(0.01, 0)), dim=c(2,2)))
  
  # #1 is dropout, #2 below threshold (imputation occurs)
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(1, 0), dim=c(1, N))
  data <- array(c(1.2, 4.5), dim=c(1,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(0, 0)), dim=c(2,2)))
  
  # #1 is dropout, #2 above threshold
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(1, 0), dim=c(1, N))
  data <- array(c(1.2, 5.2), dim=c(1,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(16, 0)), dim=c(2,2)))
  
  # 3 features
  dist <- array(0, dim=c(N, N))
  dropouts <- array(c(c(0, 0, 0),
                      c(1, 1, 0)), dim=c(3, N))
  data <- array(c(c(4.2, 5.7, 4.5),
                  c(1.2, 0.7, 2.5)), dim=c(3,N))
  threshold <- 5.0
  cpp_dist(dist, dropouts, data, N, threshold)
  expect_equal(dist, array(c(c(0, 0), c(29, 0)), dim=c(2,2)))
})
