library(tinytest)
library(repmix)

## w = x/(x + y) should give the same marginal likelihood
m1 <- marglik(tr = 0.09, sr = 0.05, to = 0.21, so = 0.05, x = 2, y = 2, m = 0, v = 2)
m2 <- marglik(tr = 0.09, sr = 0.05, to = 0.21, so = 0.05, w = 0.5, m = 0, v = 2)
expect_equal(m1, m2)

## TODO add more tests
