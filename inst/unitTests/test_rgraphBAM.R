## http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/
## checkEquals: Are two objects equal, including named attributes?
## checkEqualsNumeric: Are two numeric values equal?
## checkIdentical: Are two objects exactly the same?
## checkTrue: Does an expression evaluate to TRUE?
## checkException: Does an expression raise an error?

## rgraphBAM()
test_rgraphBAM <- function() {

  ## check that random d-regular graphs have a constant degree
  set.seed(1234)
  g <- rgraphBAM(dRegularMarkedGraphParam(pI=2, pY=10, d=3))
  d <- graph::degree(g)
  checkIdentical(unique(d), 3)
}
