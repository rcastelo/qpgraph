\name{qpgraph-deprecated}
\alias{qpgraph-deprecated}
\title{Deprecated Functions in Package \pkg{qpgraph}}
\description{
  These functions are provided for compatibility with older versions of
  \pkg{qpgraph}  only, and may be defunct as soon as the next release.
}
\usage{
qpRndHMGM(nDiscrete=1, nContinuous=3, d=2, mixedIntStrength=5, rho=0.5, G=NULL)
qpSampleFromHMGM(n=10, hmgm=qpRndHMGM())
}
\details{
  \code{qpRndHMGM()} has been replaced by \code{rHMgmm()}.

  \code{qpSampleFromHMGM()} has been replaced by \code{rcmvnorm()}.

  The original help page for these functions is available as
  \code{help("oldName-deprecated")} (note the quotes).
}
\seealso{
  \code{\link{Deprecated}}
}
\keyword{internal}
\keyword{misc}
