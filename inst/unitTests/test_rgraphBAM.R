test_rgraphBAM <- function() {
  set.seed(1234)
  g <- rgraphBAM(dRegularMarkedGraphParam(pI=2, pY=10, d=3))
  d <- graph::degree(g)
  checkIdentical(unique(d), 3)
}
