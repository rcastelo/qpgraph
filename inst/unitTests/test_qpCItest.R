## http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/
## checkEquals: Are two objects equal, including named attributes?
## checkEqualsNumeric: Are two numeric values equal?
## checkIdentical: Are two objects exactly the same?
## checkTrue: Does an expression evaluate to TRUE?
## checkException: Does an expression raise an error?

## qpCItest()
test_qpCItest <- function() {

  nObs <- 100

  ## the following adjacency matrix describes an undirected graph
  ## where vertex 3 is conditionally independent of 4 given 1 AND 2
  A <- matrix(c(FALSE,  TRUE,  TRUE,  TRUE,
                TRUE,  FALSE,  TRUE,  TRUE,
                TRUE,   TRUE, FALSE, FALSE,
                TRUE,   TRUE, FALSE, FALSE), nrow=4, ncol=4, byrow=TRUE)

  set.seed(1234)
  gmm <- rUGgmm(A)
  K <- solve(gmm$sigma)

  ## check that the missing edge corresponds to a zero in the corresponding
  ## cell of the inverse covariance matrix
  checkEqualsNumeric(K[!A & upper.tri(A)], 0,
                     msg="Missing edges match zero pattern in the inverse covariance",
                     tolerance=1e-10)
  X <- rmvnorm(nObs, gmm)
  cit <- qpCItest(X, i=3, j=4, Q=1:2, long.dim.are.variables=FALSE)

  ## check that the p-value is non-significant and about 0.44
  checkEqualsNumeric(cit$p.value, 0.44, msg="Conditional indepencence test for pure continuous data",
                     tolerance=0.01)
}
