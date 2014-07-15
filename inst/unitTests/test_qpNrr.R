## http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/
## checkEquals: Are two objects equal, including named attributes?
## checkEqualsNumeric: Are two numeric values equal?
## checkIdentical: Are two objects exactly the same?
## checkTrue: Does an expression evaluate to TRUE?
## checkException: Does an expression raise an error?

## qpNrr()
test_qpNrr <- function() {

  set.seed(123)
              
  ## simulate an undirected Gaussian graphical model
  ## with no edges
  model <- rUGgmm(dRegularGraphParam(p=50, d=0), rho=0.5)

  ## simulate data from this model
  X <- rmvnorm(n=100, model)

  ## estimate non-rejection rates with q=3
  nrr.estimates <- qpNrr(X, q=1, identicalQs=FALSE, long.dim.are.variables=FALSE, verbose=FALSE)

  ## create an adjacency matrix of the undirected graph
  ## determining the undirected Gaussian graphical model
  A <- as(model$g, "matrix") == 1

  ## mean value of the non-rejection rates for missing pure continuous edges
  checkEqualsNumeric(mean(nrr.estimates[upper.tri(nrr.estimates) & !A]), 0.95,
                     msg=" non-rejection rate with pi_ij=1 should approach alpha=0.05",
                     tolerance=0.01)

  map <- qtl::sim.map(len=100,          ## genetic map consisting of one chromosome 100 cM long
                      n.mar=10,         ## with 10 markers equally spaced along the chromosome
                      anchor.tel=FALSE,
                      eq.spacing=TRUE)

  ## simulate an eQTL network with one cis-eQTL and no gene-gene associations
  sim.eqtl <- reQTLcross(eQTLcrossParam(map=map, genes=50, cis=1L,
                                        networkParam=dRegularGraphParam(d=0)),
                                        rho=0.5, a=1.0)

  ## simulate data (a qtl/cross object) using this eQTL network
  cross <- sim.cross(map, sim.eqtl, n.ind=100)

  eQTLgene <- alleQTL(sim.eqtl)$gene

  ## calculate all pairwise NRR fixing the eQTL gene in each conditioning set
  nrr.estimates <- qpNrr(cross, q=2, restrict.Q=setdiff(sim.eqtl$model$Y, eQTLgene),
                         fix.Q=eQTLgene, verbose=FALSE)

  ## since the eQTL gene was fixed by conditioning, no NRR has been calculated for this
  ## gene, and therefore, all NRRs should approach 0.95
  checkEqualsNumeric(mean(nrr.estimates[upper.tri(nrr.estimates)], na.rm=TRUE), 0.95,
                     msg=" non-rejection rate with pi_ij=1 should approach alpha=0.05",
                     tolerance=0.01)
}
