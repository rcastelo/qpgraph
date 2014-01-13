## qpgraph package - this R code implements functions to learn qp-graphs from
##                   data, to estimate partial correlations, simulate undirected Gaussian
##                   graphical models and deal with microarray and genetic data in order
##                   to build network models of molecular regulation
##
## Copyright (C) 2013 R. Castelo and A. Roverato, with contributions from Inma Tur.
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, you can obtain one via WWW at
## http://www.gnu.org/copyleft/gpl.html, or by writing to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.



## function: qpNrr
## purpose: estimate non-rejection rates for every pair of variables
## parameters: X - data set from where to estimate the non-rejection rates
##             q - partial-correlation order to be employed
##             I - indexes or names of the variables in X that are discrete
##             restrict.Q - indexes or names of variables to which the conditioning
##                          subsets Q should be restricted. this can be a logical
##                          squared matrix indicating differerent restriction subsets
##                          per variable row-wise
##             fix.Q - indexes or names of variables that should be fixed within
##                     every conditioning subset Q
##             nTests - number of tests for each pair of variables
##             alpha - significance level of each test (Type-I error probability)
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             verbose - show progress on the calculations
##             identicalQs - use identical conditioning subsets for all pairs
##                           of variables
##             exact.test - employ an exact test when I!=NULL
##             R.code.only - flag set to FALSE when using the C implementation
##             clusterSize - size of the cluster of processors to do calculations
##                           in parallel via 'snow' and 'rlecuyer'
## return: a matrix with the estimates of the non-rejection rates

## X comes as an ExpressionSet object
setMethod("qpNrr", signature(X="ExpressionSet"),
          function(X, q=1, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, verbose=TRUE,
                   identicalQs=TRUE, exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                   nAdj2estimateTime=10) {

            use <- match.arg(use)

            p <- as.integer(nrow(X))
            h <- as.integer(ncol(Biobase::pData(X)))
            pNames <- colnames(Biobase::pData(X))

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            XP <- matrix(NA, nrow=ncol(X), ncol=0)
            I <- NULL
            if (h > 0) { ## if there are phenotypic variables, they are allowed to
                         ## to be included in pairup.i, pairup.j or fix.Q
              if (is.character(pairup.i)) {
                mt <- match(pairup.i, pNames)
                for (i in mt[!is.na(mt)]) {
                  x <- Biobase::pData(X)[, i]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[i])
                }
              } else {
                for (i in which(pairup.i > p)) {
                  x <- Biobase::pData(X)[, pairup.i[i]-p]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[pairup.i[i]-p])
                  pairup.i[i] <- p+ncol(XP)
                }
              }
              if (is.character(pairup.j)) {
                mt <- match(pairup.j, pNames)
                for (i in mt[!is.na(mt)]) {
                  x <- Biobase::pData(X)[, i]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[i])
                }
              } else {
                for (i in which(pairup.j > p)) {
                  x <- Biobase::pData(X)[, pairup.j[i]-p]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[pairup.j[i]-p])
                  pairup.j[i] <- p+ncol(XP)
                }
              }
              if (is.character(fix.Q)) {
                mt <- match(fix.Q, pNames)
                for (i in mt[!is.na(mt)]) {
                  x <- Biobase::pData(X)[, i]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[i])
                }
              } else {
                for (i in which(fix.Q > p)) {
                  x <- Biobase::pData(X)[, fix.Q[i]-p]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[fix.Q[i]-p])
                  fix.Q[i] <- p+ncol(XP)
                }
              }
            } ## end if (h > 0)

            X <- t(Biobase::exprs(X))
            X <- cbind(X, XP)
            qpgraph:::.qpNrr(X, q, I, restrict.Q, fix.Q, nTests, alpha,
                             pairup.i, pairup.j, verbose, identicalQs,
                             exact.test, use, tol, R.code.only, clusterSize,
                             startTime, nAdj2estimateTime)
          })

## X comes as a cross object
setMethod("qpNrr", signature(X="cross"),
          function(X, q=1, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, verbose=TRUE,
                   identicalQs=TRUE, exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                   nAdj2estimateTime=10) {

            use <- match.arg(use)

            p <- as.integer(qtl::nphe(X) + qtl::totmar(X))
            n <- qtl::nind(X)
            phenoNames <- colnames(X$pheno)
            markerNames <- qtl::markernames(X)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            XP <- matrix(NA_real_, nrow=n, ncol=p, dimnames=list(NULL, c(markerNames, phenoNames)))
            phclass <- sapply(X$pheno, class)
            discreteMask <- phclass == "character" | phclass == "factor" | phclass == "logical"
            I <- c(markerNames, phenoNames[discreteMask])
            XP[, markerNames] <- do.call(cbind, lapply(X$geno, function(x) x$data)) ## assuming all genotype alleles come as 1, 2, ...
            if (any(phclass == "numeric"))
              XP[, phenoNames[phclass == "numeric"]] <- as.matrix(X$pheno[, phenoNames[phclass == "numeric"]])
            if (any(phclass == "integer"))
              XP[, phenoNames[phclass == "integer"]] <- do.call("mode<-",
                                                                list(as.matrix(X$pheno[, phenoNames[phclass == "integer"]]),
                                                                     "numeric"))
            if (any(discreteMask))
              XP[, phenoNames[discreteMask]] <- as.matrix(as.data.frame(lapply(lapply(X$pheno[, phenoNames[discreteMask]], factor), as.numeric)))
            
            ## different behavior for cross objects: by default only NRR values between markers and phenotypes,
            ## and continuous phenotypes, are estimated
            if (is.null(pairup.i) && is.null(pairup.j)) {
              pairup.i <- c(markerNames, phenoNames)
              pairup.j <- setdiff(phenoNames, I)
            }

            ## different behavior for cross objects: by default the subsets Q are samples from continuous phenotypes
            if (is.null(restrict.Q))
              restrict.Q <- setdiff(phenoNames, I)

            qpgraph:::.qpNrr(XP, q, I, restrict.Q, fix.Q, nTests, alpha,
                             pairup.i, pairup.j, verbose, identicalQs,
                             exact.test, use, tol, R.code.only, clusterSize,
                             startTime, nAdj2estimateTime)
          })

## X comes as a data frame
setMethod("qpNrr", signature(X="data.frame"),
          function(X, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                   long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                   R.code.only=FALSE, clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10) {

            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(m)))
              colnames(X) <- 1:ncol(X)

            qpgraph:::.qpNrr(X, q, I, restrict.Q, fix.Q, nTests, alpha,
                             pairup.i, pairup.j, verbose, identicalQs,
                             exact.test, use, tol, R.code.only, clusterSize,
                             startTime, nAdj2estimateTime)
          })

          
## X comes as a matrix
setMethod("qpNrr", signature(X="matrix"),
          function(X, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                   long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                   R.code.only=FALSE, clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10) {

            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)

            qpgraph:::.qpNrr(X, q, I, restrict.Q, fix.Q, nTests, alpha,
                             pairup.i, pairup.j, verbose, identicalQs,
                             exact.test, use, tol, R.code.only, clusterSize,
                             startTime, nAdj2estimateTime)
          })

.qpNrr <- function(X, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, verbose=TRUE,
                   identicalQs=TRUE, exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE, clusterSize=1, startTime, nAdj2estimateTime=10) {

  if (use == "em" && !R.code.only)
    stop("use=\"em\" does not work yet with R.code.only=FALSE\n")

  cl <- NULL
 
  if (class(clusterSize)[1] == "numeric" || class(clusterSize)[1] == "integer") {
    if (clusterSize > 1) {
      ## copying ShortRead's strategy, 'get()' are to quieten R CMD check, and for no other reason
      makeCl <- get("makeCluster", mode="function")
      clSetupRNG <- get("clusterSetupRNG", mode="function")
      clEvalQ <- get("clusterEvalQ", mode="function")
      clExport <- get("clusterExport", mode="function")
      clApply <- get("clusterApply", mode="function")
      stopCl <- get("stopCluster", mode="function")
      clCall <- get("clusterCall", mode="function")
      clOpt <- get("getClusterOption", mode="function")

      if (startTime["elapsed"] == 0)
        message("Estimating non-rejection rates using a ", clOpt("type"),
                " cluster of ", clusterSize, " nodes\n")
      else
        message("Estimating time of calculation of non-rejection rates using a ", clOpt("type"),
                " cluster of ", clusterSize, " nodes\n")

      cl <- makeCl(clusterSize, type="MPI", snowlib=system.file(package="qpgraph"))
      clSetupRNG(cl)
      res <- clEvalQ(cl, require(qpgraph, quietly=TRUE))
      if (!all(unlist(res))) {
        stopCl(cl)
        stop("The package 'qpgraph' could not be loaded in some of the nodes of the cluster")
      }
      assign("clusterSize", clusterSize, envir=.GlobalEnv)
      clExport(cl, list("clusterSize"))
      rm("clusterSize", envir=.GlobalEnv)
      clApply(cl, 1:clusterSize, function(x) assign("clusterRank", x, envir=.GlobalEnv))
    }
  } else {
    if (!is.na(match("cluster", class(clusterSize))))
      cl <- clusterSize
    else
      stop("Unknown class for argument clusterSize:", class(clusterSize))
  }

  ## X the matrix of data with columns as r.v. and rows as multivariate observations

  var.names <- colnames(X)
  n.var <- ncol(X)
  N <- nrow(X)

  ## check that the q, nTests and alpha parameters have proper values

  if (q > n.var - 2)
    stop(paste("q=",q," > p-2=", n.var-2))

  if (q < 0)
    stop(paste("q=",q," < 0"))

  if (q > N - 3)
    stop(paste("q=",q," > n-3=",N-3))

  nTests <- as.integer(nTests)
  if (nTests < 1)
    stop(paste("nTests=",nTests," < 1"))

  if (alpha < 0.0 || alpha > 1.0) {
    stop(sprintf("significance level alpha is %.2f and it should lie in the interval [0,1]\n",alpha))
  }

  ## check whether there are discrete variables and whether they're properly set

  Y <- NULL
  if (!is.null(I)) {
    if (is.character(I)) {
      if (any(is.na(match(I, var.names))))
        stop("Some variables in I do not form part of the variable names of the data in X\n")
      I <- match(I, var.names)
    } else {
      if (any(is.na(match(I, 1:n.var))))
        stop("Some variables in I do not form part of the variables of the data in X\n")
    }
    Y <- (1:n.var)[-I]
  }

  ## should the following error messages stop the cluster if it has been started ??

  if (!is.null(restrict.Q)) {
    if (!is.matrix(restrict.Q) && !is.integer(restrict.Q) && !is.character(restrict.Q))
      stop("restrict.Q should be either a matrix or a vector of indexes or variables names\n")

    if (is.matrix(restrict.Q)) {
      if (is.null(I))
        stop("restrict.Q as a matrix can only be employed for restricting conditioning of discrete variables\n")

      if (rownames(restrict.Q) != colnames(X) || colnames(restrict.Q) != colnames(X))
        stop("row and column names in restrict.Q should coincide to the column names in X\n")

      if (q > min(rowSums(restrict.Q))-2)
        stop("The minimum number of variables in restrict.Q from where subsets Q of size q should be sampled is < (q+2)\n")
    } else {
      if (is.character(restrict.Q)) {
        if (any(is.na(match(restrict.Q, var.names))))
          stop("Some variables in restrict.Q do not form part of the variable names of the data in X\n")
        restrict.Q <- match(restrict.Q, var.names)
      } else {
        if (any(is.na(match(restrict.Q, 1:n.var))))
          stop("Some variables in restrict.Q do not form part of the variables of the data in X\n")
      }

      if (q > length(restrict.Q)-2)
        stop("The number of variables in restrict.Q from where subsets Q of size q should be sampled is < (q+2)\n")
    }
  }

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (is.null(pairup.i))
    pairup.i <- 1:n.var
  else {
    if (is.character(pairup.i)) {
      if (any(is.na(match(pairup.i, var.names))))
        stop("Some variables in pairup.i do not form part of the variable names of the data in X\n")
      pairup.i <- match(pairup.i, var.names)
    }
  }

  if (is.null(pairup.j)) {
    pairup.j <- 1:n.var
    if (!is.null(I)) { ## by now, interactions between discrete variables are not considered
        pairup.j <- (1:n.var)[-I]
    }
  } else {
    if (is.character(pairup.j)) {
      if (any(is.na(match(pairup.j, var.names))))
        stop("Some variables in pairup.j do not form part of the variable names of the data in X\n")
      pairup.j <- match(pairup.j, var.names)
    }
  }

  if (!is.null(fix.Q)) {
    if (q <= length(fix.Q))
      stop("q should be larger than the number of variables in fix.Q\n")

    if (is.character(fix.Q)) {
      if (any(is.na(match(fix.Q, var.names))))
        stop("Some variables in fix.Q do not form part of the variable names of the data\n")
      fix.Q <- match(fix.Q, var.names)
    } else {
      if (any(is.na(match(fix.Q, 1:n.var))))
        stop("Some variables in fix.Q do not form part of the variables of the data\n")
    }

    if (is.null(restrict.Q))
      restrict.Q <- setdiff(1:n.var, fix.Q)
    else {
      if (is.matrix(restrict.Q)) {
        if (any(apply(restrict.Q[-fix.Q, ], 1, function(x, y) intersect(x, y)) > 0))
          stop("The subsets restrict.Q and fix.Q should be disjoint.\n")
      } else {
        if (length(intersect(restrict.Q, fix.Q)) > 0)
          stop("The subsets restrict.Q and fix.Q should be disjoint.\n")
      }
    }

    ## variables in fix.Q are removed from the pairs for which nrr values are estimated
    pairup.i <- setdiff(pairup.i, fix.Q)
    pairup.j <- setdiff(pairup.j, fix.Q)
  }

  ## pair the two sets pairup.i and pairup.j without pairing the same variable
  l.pairup.i <- length(pairup.i)
  l.pairup.j <- length(pairup.j)
  l.int <- length(intersect(pairup.i, pairup.j))
  l.pairup.i.noint <- l.pairup.i - l.int
  l.pairup.j.noint <- l.pairup.j - l.int
  n.adj <- l.int * l.pairup.j.noint + l.int * l.pairup.i.noint +
           l.pairup.i.noint * l.pairup.j.noint + l.int * (l.int - 1) / 2

  pairup.ij.int <- intersect(pairup.i, pairup.j)
  pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
  pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

  nrrMatrix <- NULL

  ## estimate the actual number of necessary tests for number required by the user
  if (identicalQs && is.null(I)) {
    fractionValidQs <- 1-phyper(0, 2, n.var-2-length(fix.Q), q-length(fix.Q), lower.tail=FALSE)
    if (fractionValidQs < 0.9) {
      warning(paste(sprintf("With p=%d and q=%d the estimated fraction of valid Q sets is %.2f.", n.var, q, fractionValidQs),
                      "Increasing nTests from", nTests, "to", floor(nTests/fractionValidQs), "in order to achieve the desired precision\n", sep=" "))
      nTests <- floor(nTests / fractionValidQs)
    }
  }

  if (!R.code.only) {
    elapsedTime <- 0
    if (startTime["elapsed"] > 0) {
      elapsedTime <- (proc.time() - startTime)["elapsed"]
      startTime <- proc.time()
    }

    if (is.null(cl)) { ## single-processor execution

      if (identicalQs && is.null(I))
        nrrMatrix <- qpgraph:::.qpFastNrrIdenticalQs(X, q, restrict.Q, fix.Q,
                                                     nTests, alpha, pairup.i.noint,
                                                     pairup.j.noint, pairup.ij.int, verbose,
                                                     startTime["elapsed"], nAdj2estimateTime)
      else
        nrrMatrix <- qpgraph:::.qpFastNrr(X, I, Y, q, restrict.Q, fix.Q, nTests, alpha,
                                          pairup.i.noint, pairup.j.noint,
                                          pairup.ij.int, exact.test, verbose,
                                          startTime["elapsed"], nAdj2estimateTime)

      if (startTime["elapsed"] == 0)
        nrrMatrix <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                         Dimnames=list(var.names, var.names), x=nrrMatrix)

    } else {           ## use a cluster !
      clCall <- get("clusterCall", mode="function")
      nrrIdx <- list()
      if (verbose && startTime["elapsed"] == 0) { ## no cluster progress-call when only estimating time
        if (identicalQs && is.null(I))
          nrrIdx <- clPrCall(cl, qpgraph:::.qpFastNrrIdenticalQsPar, n.adj, X,
                             q, restrict.Q, fix.Q, nTests, alpha, pairup.i.noint,
                             pairup.j.noint, pairup.ij.int, verbose, FALSE,
                             nAdj2estimateTime)
        else
          nrrIdx <- clPrCall(cl, qpgraph:::.qpFastNrrPar, n.adj, X, I, Y, q,
                             restrict.Q, fix.Q, nTests, alpha, pairup.i.noint,
                             pairup.j.noint, pairup.ij.int, exact.test, verbose,
                             FALSE, nAdj2estimateTime)
      } else {
        if (identicalQs && is.null(I))
          nrrIdx <- clCall(cl, qpgraph:::.qpFastNrrIdenticalQsPar, X, q,
                           restrict.Q, fix.Q, nTests, alpha, pairup.i.noint,
                           pairup.j.noint, pairup.ij.int, verbose,
                           startTime["elapsed"] > 0, nAdj2estimateTime)
        else
          nrrIdx <- clCall(cl, qpgraph:::.qpFastNrrPar, X, I, Y, q, restrict.Q,
                           fix.Q, nTests, alpha, pairup.i.noint, pairup.j.noint,
                           pairup.ij.int, exact.test, verbose,
                           startTime["elapsed"] > 0, nAdj2estimateTime)
      }

      if (startTime["elapsed"] > 0) {
        ## the following calculation makes important part of the estimation of the time
        ## it assumes that the estimated time per processor is stored on the first position of 'nrr'
        ## and uses the median of the times estimated for each processor to try to be robust against
        ## fluctuations on the execution time taken in some processors
        elapsedTime <- elapsedTime + median(sapply(nrrIdx, function(x) x$nrr[1]))
        startTime <- proc.time()
      }

      if (class(clusterSize)[1] == "numeric" || class(clusterSize)[1] == "integer")
        stopCl(cl)

      nrrMatrix <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                       Dimnames=list(var.names, var.names),
                       x=rep(as.double(NA), n.var*(n.var-1)/2+n.var)) 
      nrrMatrix@x[do.call("c", lapply(nrrIdx, function(x) x$idx))] <-
        do.call("c", lapply(nrrIdx, function(x) x$nrr))

      if (startTime["elapsed"] > 0) {
        elapsedTime <- elapsedTime + (proc.time() - startTime)["elapsed"]
        d <- as.vector(floor(elapsedTime / (24*3600)))
        h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
        m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
        s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
        nrrMatrix <- c(days=d, hours=h, minutes=m, seconds=s)
      }
    }

    return(nrrMatrix)
  }

  if (identicalQs && is.null(I)) {
    nrrMatrix <- .qpNrrIdenticalQs(X, q, restrict.Q, fix.Q, nTests, alpha,
                                   pairup.i.noint, pairup.j.noint, pairup.ij.int,
                                   verbose, startTime, nAdj2estimateTime)

    return(nrrMatrix)
  }

  S <- ssd <- mapX2ssd <- NULL

  missingMask <- apply(X, 1, function(x) any(is.na(x)))
  missingData <- any(missingMask)

  if (!missingData) {
    if (!is.null(I)) {  ## calculate the uncorrected sum of squares and deviations
      ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
      mapX2ssd <- match(var.names, colnames(ssd))
      ## names(mapX2ssd) <- colnames(X) ## is this necessary
    } else             ## calculate sample covariance matrix
      S <- qpCov(X)
  }

  if (!is.null(I)) {
    nLevels <- rep(NA_integer_, times=ncol(X))
    nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
    if (any(nLevels[I] == 1))
      stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))
  }

  ## the idea is to return an efficiently stored symmetric matrix
  nrrMatrix <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                   Dimnames=list(var.names, var.names),
                   x=rep(as.double(NA), n.var*(n.var-1)/2+n.var))

  elapsedTime <- 0
  if (startTime["elapsed"] > 0) {
    elapsedTime <- (proc.time() - startTime)["elapsed"]
    startTime <- proc.time()
  }

  ppct <- -1
  k <- 0
  pb <- NULL
  if (verbose && elapsedTime == 0)
    pb <- txtProgressBar(style=3)
  rQs <- NULL
  if (!is.null(restrict.Q) && !is.matrix(restrict.Q))
    rQs <- restrict.Q

  nrr <- NA

  ## intersection variables against ij-exclusive variables
  for (i in pairup.ij.int) {
    for (j in c(pairup.i.noint,pairup.j.noint)) {

      if (is.null(I))
        nrr <- qpgraph:::.qpEdgeNrr(X, S, i, j, q, rQs, fix.Q, nTests,
                                    alpha, R.code.only=TRUE)
      else {
        if (!is.null(restrict.Q) && is.matrix(restrict.Q))
            rQs <- union(which(restrict.Q[i, ]), which(restrict.Q[j, ]))

        nrr <- qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q, rQs, fix.Q,
                                        nTests, alpha, exact.test, use, tol, R.code.only=TRUE)
      }

      nrrMatrix[j,i] <- nrrMatrix[i,j] <- nrr
      k <- k + 1
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose && elapsedTime == 0) {
        setTxtProgressBar(pb, pct/100)
        ppct <- pct
      }
    }
    if (elapsedTime > 0 && k == nAdj2estimateTime)
      break;
  }

  ## i-exclusive variables against j-exclusive variables
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in pairup.i.noint) {
      for (j in pairup.j.noint) {

        if (is.null(I))
          nrr <- qpgraph:::.qpEdgeNrr(X, S, i, j, q, rQs, fix.Q, nTests,
                                      alpha, R.code.only=TRUE)
        else {
          if (!is.null(restrict.Q) && is.matrix(restrict.Q))
            rQs <- union(which(restrict.Q[i, ]), which(restrict.Q[j, ]))

          nrr <- qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q, rQs, fix.Q,
                                          nTests, alpha, exact.test, use, tol, R.code.only=TRUE)
        }

        nrrMatrix[j,i] <- nrrMatrix[i,j] <- nrr
        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  ## intersection variables against themselves (avoiding pairing of the same)
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in seq(along=pairup.ij.int[-1])) { ## 1:(l.int-1) does not work if l.int == 0
      i2 <- pairup.ij.int[i]

      for (j in (i+1):l.int) {
        j2 <- pairup.ij.int[j]

        if (is.null(I))
          nrr <- qpgraph:::.qpEdgeNrr(X, S, i2, j2, q, rQs, fix.Q, nTests,
                                      alpha, R.code.only=TRUE)
        else {
          if (!is.null(restrict.Q) && is.matrix(restrict.Q))
            rQs <- union(which(restrict.Q[i2, ]), which(restrict.Q[j2, ]))

          nrr <- qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i2, j2, q, rQs, fix.Q,
                                          nTests, alpha, exact.test, use, tol, R.code.only=TRUE)

        }

        nrrMatrix[j2,i2] <- nrrMatrix[i2,j2] <- nrr
        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  if (verbose && elapsedTime == 0) {
    close(pb)
  }

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + ((proc.time()-startTime)["elapsed"]/k) * n.adj
    startTime <- proc.time()
  }

  ## this is necessary till we find out how to properly assign values in a dspMatrix
  nrrMatrix <- as(nrrMatrix, "dspMatrix")

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + (proc.time()-startTime)["elapsed"]
    d <- as.vector(floor(elapsedTime / (24*3600)))
    h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
    m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
    s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
    nrrMatrix <- c(days=d, hours=h, minutes=m, seconds=s)
  }

  return(nrrMatrix)
}

.qpNrrIdenticalQs <- function(X, q, restrict.Q, fix.Q, nTests, alpha,
                              pairup.i.noint, pairup.j.noint, pairup.ij.int,
                              verbose, startTime, nAdj2estimateTime) {

  ## X the matrix of data with columns as r.v. and rows as multivariate observations
  var.names <- colnames(X)
  n.var <- ncol(X)
  N <- nrow(X)

  ## how many adjacencies do we have to calculate
  l.int <- length(pairup.ij.int)
  l.pairup.i.noint <- length(pairup.i.noint)
  l.pairup.j.noint <- length(pairup.j.noint)
  n.adj <- l.int * l.pairup.j.noint + l.int * l.pairup.i.noint +
           l.pairup.i.noint * l.pairup.j.noint + l.int * (l.int - 1) / 2

  ## calculate sample covariance matrix
  S <- qpCov(X)

  ## sample the Q sets and pre-calculate the inverse matrices
  if (is.null(restrict.Q))
    restrict.Q <- 1:n.var
  n.fQ <- length(fix.Q)
  restrict.Q <- setdiff(restrict.Q, fix.Q)

  ## TODO: work with complete.obs on the margins
  Qs <- as.list(array(dim=nTests))
  Qs <- lapply(Qs, function(x, rQ, fQ, n.fQ) c(sample(rQ, size=q-n.fQ, replace=FALSE), fQ),
               restrict.Q, fix.Q, n.fQ)
  S22invs <- lapply(Qs, function(x) solve(S[x, x]) )

  ## the idea is to return an efficiently stored symmetric matrix
  nrrMatrix <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                   Dimnames=list(var.names, var.names),
                   x=rep(as.double(NA), n.var*(n.var-1)/2+n.var))

  elapsedTime <- 0
  if (startTime["elapsed"] > 0) {
    elapsedTime <- (proc.time() - startTime)["elapsed"]
    startTime <- proc.time()
  }

  ppct <- -1
  k <- 0
  pb <- NULL
  if (verbose && elapsedTime == 0)
    pb <- txtProgressBar(style=3)

  ## intersection variables against ij-exclusive variables
  for (i in pairup.ij.int) {
    for (j in c(pairup.i.noint,pairup.j.noint)) {
      nrrMatrix[j,i] <- nrrMatrix[i,j] <-
        qpgraph:::.qpEdgeNrrIdenticalQs(S, Qs, S22invs, i, j, q, nTests,
                                        alpha, R.code.only=TRUE)
      k <- k + 1
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose && elapsedTime == 0) {
        setTxtProgressBar(pb, pct/100)
        ppct <- pct
      }
    }
    if (elapsedTime > 0 && k == nAdj2estimateTime)
      break;
  }

  ## i-exclusive variables against j-exclusive variables
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in pairup.i.noint) {
      for (j in pairup.j.noint) {
        nrrMatrix[j,i] <- nrrMatrix[i,j] <-
          qpgraph:::.qpEdgeNrrIdenticalQs(S, Qs, S22invs, i, j, q, nTests,
                                          alpha, R.code.only=TRUE)
        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  l.int <- length(pairup.ij.int)

  ## intersection variables against themselves (avoiding pairing of the same)
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in seq(along=pairup.ij.int[-1])) { ## 1:(l.int-1) does not work if l.int==0
      i2 <- pairup.ij.int[i]

      for (j in (i+1):l.int) {
        j2 <- pairup.ij.int[j]
        nrrMatrix[j2,i2] <- nrrMatrix[i2,j2] <-
          qpgraph:::.qpEdgeNrrIdenticalQs(S, Qs, S22invs, i2, j2, q, nTests,
                                          alpha, R.code.only=TRUE)
        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  if (verbose && elapsedTime == 0) {
    close(pb) 
  }

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + ((proc.time()-startTime)["elapsed"]/k) * n.adj
    startTime <- proc.time()
  }

  ## this is necessary till we find out how to properly assign values in a dspMatrix
  nrrMatrix <- as(nrrMatrix, "dspMatrix")

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + (proc.time()-startTime)["elapsed"]
    d <- as.vector(floor(elapsedTime / (24*3600)))
    h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
    m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
    s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
    nrrMatrix <- c(days=d, hours=h, minutes=m, seconds=s)
  }

  return(nrrMatrix)
}



## function: qpAvgNrr
## purpose: estimate average non-rejection rates for every pair of variables
## parameters: X - data set from where to estimate the average non-rejection
##                 rates
##             qOrders - either a number of partial-correlation orders or a
##                       vector of particular orders to be employed in the
##                       calculation
##             I - indexes or names of the variables in X that are discrete
##             restrict.Q - indexes or names of variables to which the conditioning
##                          subsets Q should be restricted. this can be a logical
##                          squared matrix indicating differerent restriction subsets
##                          per variable row-wise
##             fix.Q - indexes or names of variables that should be fixed within
##                     every conditioning subset Q
##             nTests - number of tests to perform for each pair of variables
##             alpha - significance level of each test (Type-I error probability)
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             type - type of average (by now only the arithmetic mean is
##                    available)
##             verbose - show progress of the calculations
##             identicalQs - use identical conditioning subsets for all pairs
##                           of variables
##             exact.test - employ an exact test when I!=NULL
##             R.code.only - flag set to FALSE when using the C implementation
##             clusterSize - size of the cluster of processors to do calculations
##                           in parallel via 'snow' and 'rlecuyer'
## return: a matrix with the estimates of the average non-rejection rates

setGeneric("qpAvgNrr", function(X, ...) standardGeneric("qpAvgNrr"))

## X comes as an ExpressionSet object
setMethod("qpAvgNrr", signature(X="ExpressionSet"),
          function(X, qOrders=4, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, type=c("arith.mean"),
                   verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                   clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10) {

            type <- match.arg(type)
            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            X_I <- NULL
            if (!is.null(I)) {
              if (!is.character(I))
                stop("When X is an ExpressionSet, I can only contain variable names from the associated phenotypic data.")
              if (any(is.na(match(I, Biobase::varLabels(X)))))
                stop(sprintf("%s do(es) not form part of the phenotypic data in the ExpressionSet object X.",
                             I[is.na(match(I, Biobase::varLabels(X)))]))

              X_I <- apply(Biobase::pData(X)[, I, drop=FALSE], 2,
                           function(x) as.double(as.factor(as.character(x))))
            }

            X <- cbind(t(Biobase::exprs(X)), X_I)
            qpgraph:::.qpAvgNrr(X, qOrders, I, restrict.Q, fix.Q, nTests, alpha, pairup.i,
                                pairup.j, type, verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime, nAdj2estimateTime)
          })

## X comes as a data frame
setMethod("qpAvgNrr", signature(X="data.frame"),
          function(X, qOrders=4, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, long.dim.are.variables=TRUE,
                   type=c("arith.mean"), verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE, clusterSize=1,
                   estimateTime=FALSE, nAdj2estimateTime=10) {

            type <- match.arg(type)
            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X), decreasing=TRUE, index.return=TRUE)$ix[1] == 1)
              X <- t(X)
            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpAvgNrr(X, qOrders, I, restrict.Q, fix.Q, nTests, alpha, pairup.i,
                                pairup.j, type, verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime, nAdj2estimateTime)
          })
          
## X comes as a matrix
setMethod("qpAvgNrr", signature(X="matrix"),
          function(X, qOrders=4, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                   alpha=0.05, pairup.i=NULL, pairup.j=NULL, long.dim.are.variables=TRUE,
                   type=c("arith.mean"), verbose=TRUE, identicalQs=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                   R.code.only=FALSE, clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10) {

            type <- match.arg(type)
            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            if (long.dim.are.variables &&
                sort(dim(X), decreasing=TRUE, index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpAvgNrr(X, qOrders, I, restrict.Q, fix.Q, nTests, alpha, pairup.i,
                                pairup.j, type, verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime,
                                nAdj2estimateTime)
          })

.qpAvgNrr <- function(X, qOrders=4, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                      nTests=100, alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                      type=c("arith.mean"), verbose=TRUE, identicalQs=TRUE,
                      exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                      R.code.only=FALSE, clusterSize=1, startTime, nAdj2estimateTime) {

  cl <- 1
 
  if (clusterSize > 1) {
    ## copying ShortRead's strategy, 'get()' are to quieten R CMD check, and for no other reason
    makeCl <- get("makeCluster", mode="function")
    clSetupRNG <- get("clusterSetupRNG", mode="function")
    clEvalQ <- get("clusterEvalQ", mode="function")
    clExport <- get("clusterExport", mode="function")
    clApply <- get("clusterApply", mode="function")
    stopCl <- get("stopCluster", mode="function")
    clCall <- get("clusterCall", mode="function")
    clOpt <- get("getClusterOption", mode="function")

    if (startTime["elapsed"] == 0)
      message("Estimating average non-rejection rates using a ", clOpt("type"),
              " cluster of ", clusterSize, " nodes\n")
    else
      message("Estimating time of calculation of average non-rejection rates using a ",
              clOpt("type"), " cluster of ", clusterSize, " nodes\n")

    cl <- makeCl(clusterSize, type="MPI", snowlib=system.file(package="qpgraph"))
    clSetupRNG(cl)
    res <- clEvalQ(cl, require(qpgraph, quietly=TRUE))
    if (!all(unlist(res))) {
      stopCl(cl)
      stop("The package 'qpgraph' could not be loaded in some of the nodes of the cluster")
    }
    assign("clusterSize", clusterSize, envir=.GlobalEnv)
    clExport(cl, list("clusterSize"))
    rm("clusterSize", envir=.GlobalEnv)
    clApply(cl, 1:clusterSize, function(x) assign("clusterRank", x, envir=.GlobalEnv))
    clApply(cl, 1:clusterSize, function(x) assign("clusterRank", x, envir=.GlobalEnv))
  }

  var.names <- colnames(X)
  n.var <- ncol(X)
  N <- nrow(X)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (length(qOrders) == 1) {
    if (qOrders > min(n.var, N) - 3)
      stop(sprintf("qOrders=%d is larger than the number of available q-orders for the given data set (%d)\n",
                   qOrders, min(n.var, N) - 3))

    qOrders <- as.integer(round(seq(1, min(n.var, N) - 3,
                                    by=(min(n.var, N) - 3) / qOrders), digits=0))
  } else {
    qOrders <- as.integer(qOrders)
    if (min(qOrders) < 1 || max(qOrders) > min(n.var-3, N-3))
      stop(sprintf("for the given data set q-orders should lie in the range [1,%d]\n",
                   min(n.var-3, N-3)))
  }

  w <- 1 / length(qOrders)
  ## the idea is to return an efficiently stored symmetric matrix
  avgNrrMatrix <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                      Dimnames=list(var.names, var.names),
                      x=rep(as.double(0), n.var*(n.var-1)/2+n.var))

  elapsedTime <- 0

  for (q in qOrders) {
    if (verbose && startTime["elapsed"] == 0)
      cat(sprintf("q=%d\n",q))

    thisNrr <- qpgraph:::.qpNrr(X, q, I, restrict.Q, fix.Q, nTests, alpha, pairup.i,
                                pairup.j, verbose, identicalQs, exact.test, use, tol,
                                R.code.only, cl, startTime, nAdj2estimateTime)

    if (startTime["elapsed"] > 0) {
      elapsedTime <- elapsedTime + thisNrr["days"]*24*3600 + thisNrr["hours"]*3600 +
                     thisNrr["minutes"]*60 + thisNrr["seconds"]
      startTime <- proc.time()
    } else {
      ## this is necessary till we find out how to sum two dspMatrices getting a dspMatrix
      avgNrrMatrix <- as(avgNrrMatrix + w * thisNrr, "dspMatrix")
    }
  }

  if (clusterSize > 1 && !is.null(cl)) {
    stopCl(cl)

    elapsedTime <- elapsedTime + (proc.time() - startTime)["elapsed"]
  }

  if (startTime["elapsed"] > 0) {
    d <- as.vector(floor(elapsedTime / (24*3600)))
    h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
    m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
    s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
    avgNrrMatrix <- c(days=d, hours=h, minutes=m, seconds=s)
  }

  return(avgNrrMatrix)
}



## function: qpGenNrr
## purpose: estimate average non-rejection rates for every pair of variables
## parameters: X - data set from where to estimate the average non-rejection
##                 rates
##             datasetIdx - index vector of the different datasets. if it is
##                          a single number, it indicates the column in the
##                          phenotypic data (ExpressionSet) or in the data
##                          frame or in the matrix that indicates what samples
##                          belong to what dataset. if it is a name then it
##                          indicates the names of the phenotipic variable with
##                          this information. if it is a vector with as many
##                          positions as samples, then it contains itself the
##                          indexes about what sample belongs to what dataset
##             qOrders - either NULL indicating that a default guess on the q
##                       order will be performed for each data set or a
##                       vector of particular orders to be employed in the
##                       calculation
##             I - indexes or names of the variables in X that are discrete
##             restrict.Q - indexes or names of variables to which the conditioning
##                          subsets Q should be restricted. this can be a logical
##                          squared matrix indicating differerent restriction subsets
##                          per variable row-wise
##             fix.Q - indexes or names of variables that should be fixed within
##                     every conditioning subset Q
##             return.all - logical; set to TRUE if all intervining non-rejection
##                          rates should be return in a list; FALSE (default) if
##                          only generalized non-rejection rates should be returned
##             nTests - number of tests to perform for each pair of variables
##             alpha - significance level of each test (Type-I error probability)
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             verbose - show progress of the calculations
##             identicalQs - use identical conditioning subsets for all pairs
##                           of variables
##             exact.test - employ an exact test when I!=NULL
##             R.code.only - flag set to FALSE when using the C implementation
##             clusterSize - size of the cluster of processors to do calculations
##                           in parallel via 'snow' and 'rlecuyer'
## return: a matrix with the estimates of the average non-rejection rates

setGeneric("qpGenNrr", function(X, ...) standardGeneric("qpGenNrr"))

## maybe is better to force datasetIdx to be integer in order to get an order of the
## datasets that matches the order of the q-orders provided

## X comes as an ExpressionSet object
setMethod("qpGenNrr", signature(X="ExpressionSet"),
            function(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                     return.all=FALSE, nTests=100, alpha=0.05, pairup.i=NULL,
                     pairup.j=NULL, verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                     use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                     clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10) {

            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            if ((is.null(Biobase::pData(X)) || ncol(Biobase::pData(X)) < 1) && length(datasetIdx) != dim(X)[2])
              stop("Either supply a vector indexing the data sets to which each sample belongs to, or add a column with this information to the phenotypic data of the ExpressionSet indicating then which one is that column\n")

            if (length(datasetIdx) != dim(X)[2] && length(datasetIdx) != 1)
              stop("Argument 'datasetIdx' should be either a single number, or a character string, indicating the column of the phenotypic data of the ExpressionSet containing the indexes to the datasets. Alternatively, it can be a vector of these indexes with as many positions as samples\n")

            if (length(datasetIdx) == 1) {
              if (is.character(datasetIdx))
                datasetIdx <- match(datasetIdx, colnames(Biobase::pData(X)))
              else {
                if (is.integer(datasetIdx) || is.numeric(datasetIdx))
                  datasetIdx <- match(datasetIdx, 1:ncol(Biobase::pData(X)))
              }

              if (is.na(datasetIdx) || (!is.character(datasetIdx) && !is.integer(datasetIdx)))
                stop("Argument 'datasetIdx' does not match any phenotypic column in the input ExpressionSet X. Please look at Biobase::pData(X)\n")
            }

            if (length(datasetIdx) != dim(X)[2])
              datasetIdx <- Biobase::pData(X)[, datasetIdx]

            if (!is.null(qOrders) && any(is.na(qOrders[unique(datasetIdx)])))
              stop("Some values in 'datasetIdx' do not match any position in 'qOrders'\n")

            if (!is.null(qOrders) && is.null(names(qOrders)))
              stop("When they are specified, values in 'qOrders' should have names matching the data sets index names\n")

            X_I <- NULL
            if (!is.null(I)) {
              if (!is.character(I))
                stop("When X is an ExpressionSet, I can only contain variable names from the associated phenotypic data.")
              if (any(is.na(match(I, Biobase::varLabels(X)))))
                stop(sprintf("%s do(es) not form part of the phenotypic data in the ExpressionSet object X.",
                             I[is.na(match(I, Biobase::varLabels(X)))]))

              X_I <- apply(Biobase::pData(X)[, I, drop=FALSE], 2,
                           function(x) as.double(as.factor(as.character(x))))
            }

            X <- cbind(t(Biobase::exprs(X)), X_I)

            qpgraph:::.qpGenNrr(X, datasetIdx, qOrders, I, restrict.Q, fix.Q,
                                return.all, nTests, alpha, pairup.i, pairup.j,
                                verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime, nAdj2estimateTime)
          })

## X comes as a data frame
setMethod("qpGenNrr", signature(X="data.frame"),
          function(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                   return.all=FALSE, nTests=100, alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                   long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE, clusterSize=1,
                   estimateTime=FALSE, nAdj2estimateTime=10) {

            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)

            if (length(datasetIdx) != dim(X)[1] && length(datasetIdx) != 1)
              stop("Argument 'datasetIdx' should be either a single number, or a character string, indicating the column (or row) of the input data frame X containing the indexes to the datasets. Alternatively, it can be a vector of these indexes with as many positions as samples\n")

            if (length(datasetIdx) == 1) {
              if (is.character(datasetIdx))
                datasetIdx <- match(datasetIdx, colnames(X))
              else {
                if (is.integer(datasetIdx) || is.numeric(datasetIdx))
                  datasetIdx <- match(datasetIdx, 1:ncol(X))
              }

              if (is.na(datasetIdx) || (!is.character(datasetIdx) && !is.integer(datasetIdx) &&
                                        !is.numeric(datasetIdx)))
                stop("Argument 'datasetIdx' does not match any column (or row) in the input data frame X\n")
            }

            if (length(datasetIdx) != dim(X)[1]) {
              tmp <- X[, datasetIdx]
              X <- X[, -datasetIdx]
              datasetIdx <- tmp
            }

            if (!is.null(qOrders) && any(is.na(qOrders[unique(datasetIdx)])))
              stop("Some values in 'datasetIdx' do not match any position in 'qOrders'\n")

            if (!is.null(qOrders) && is.null(names(qOrders)))
              stop("When they are specified, values in 'qOrders' should have names matching the data sets index names\n")

            qpgraph:::.qpGenNrr(X, datasetIdx, qOrders, I, restrict.Q, fix.Q,
                                return.all, nTests, alpha, pairup.i, pairup.j,
                                verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime, nAdj2estimateTime)
          })

          
## X comes as a matrix
setMethod("qpGenNrr", signature(X="matrix"),
          function(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                   return.all=FALSE, nTests=100, alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                   long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE, clusterSize=1,
                   estimateTime=FALSE, nAdj2estimateTime=10) {

            use <- match.arg(use)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
              stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

            if (long.dim.are.variables &&
              sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)

            if (length(datasetIdx) != dim(X)[1] && length(datasetIdx) != 1)
              stop("Argument 'datasetIdx' should be either a single number, or a character string, indicating the column (or row) of the input matrix X containing the indexes to the datasets. Alternatively, it can be a vector of these indexes with as many positions as samples\n")

            if (length(datasetIdx) == 1) {
              if (is.character(datasetIdx))
                datasetIdx <- match(datasetIdx, colnames(X))
              else {
                if (is.integer(datasetIdx) || is.numeric(datasetIdx))
                  datasetIdx <- match(datasetIdx, 1:ncol(X))
              }

              if (is.na(datasetIdx) || (!is.character(datasetIdx) && !is.integer(datasetIdx) &&
                                        !is.numeric(datasetIdx)))
                stop("Argument 'datasetIdx' does not match any column (or row) in the input matrix X\n")
            }

            if (length(datasetIdx) != dim(X)[1]) {
              tmp <- X[, datasetIdx]
              X <- X[, -datasetIdx]
              datasetIdx <- tmp
            }

            if (!is.null(qOrders) && any(is.na(qOrders[unique(datasetIdx)])))
              stop("Some values in 'datasetIdx' do not match any position in 'qOrders'\n")

            if (!is.null(qOrders) && is.null(names(qOrders)))
              stop("When they are specified, values in 'qOrders' should have names matching the data sets index names\n")

            qpgraph:::.qpGenNrr(X, datasetIdx, qOrders, I, restrict.Q, fix.Q,
                                return.all, nTests, alpha, pairup.i, pairup.j,
                                verbose, identicalQs, exact.test, use, tol,
                                R.code.only, clusterSize, startTime, nAdj2estimateTime)
          })


.qpGenNrr <- function(X, datasetIdx, qOrders=NULL, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                      return.all=FALSE, nTests=100, alpha=0.05, pairup.i=NULL,
                      pairup.j=NULL, verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                      use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                      clusterSize=1, startTime, nAdj2estimateTime) {

  cl <- 1
 
  if (clusterSize > 1) {
    ## copying ShortRead's strategy, 'get()' are to quieten R CMD check, and for no other reason
    makeCl <- get("makeCluster", mode="function")
    clSetupRNG <- get("clusterSetupRNG", mode="function")
    clEvalQ <- get("clusterEvalQ", mode="function")
    clExport <- get("clusterExport", mode="function")
    clApply <- get("clusterApply", mode="function")
    stopCl <- get("stopCluster", mode="function")
    clCall <- get("clusterCall", mode="function")
    clOpt <- get("getClusterOption", mode="function")

    if (startTime["elapsed"] == 0)
      message("Estimating generalized non-rejection rates using a ", clOpt("type"),
              " cluster of ", clusterSize, " nodes\n")
    else
      message("Estimating time of calculation of generalized non-rejection rates using a ",
              clOpt("type"), " cluster of ", clusterSize, " nodes\n")

    cl <- makeCl(clusterSize, type="MPI", snowlib=system.file(package="qpgraph"))
    clSetupRNG(cl)
    res <- clEvalQ(cl, require(qpgraph, quietly=TRUE))
    if (!all(unlist(res))) {
      stopCl(cl)
      stop("The package 'qpgraph' could not be loaded in some of the nodes of the cluster")
    }
    assign("clusterSize", clusterSize, envir=.GlobalEnv)
    clExport(cl, list("clusterSize"))
    rm("clusterSize", envir=.GlobalEnv)
    clApply(cl, 1:clusterSize, function(x) assign("clusterRank", x, envir=.GlobalEnv))
  }

  var.names <- colnames(X)
  n.var <- ncol(X)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  datasetIdx <- as.character(datasetIdx)
  N <- table(datasetIdx)

  if (any(N < 3))
    stop("Dataset(s) ", paste(names(N)[which(N < 3)], collapse=","), " has/have less than 3 samples\n")

  ## when qOrders is NULL the median of the possible q-orders is taken for each dataset
  if (is.null(qOrders))
    qOrders <- floor(sapply(N, function(x) median(seq(1, x-3))))

  ## validate qOrders
  if (min(qOrders) < 1 || any(qOrders > N[names(qOrders)]-3))
    stop("Within each data set its q-order should lie in the range [1, N-3] with N being the corresponding number of samples\n")

  ## contribution of each dataset is proportional to its sample size
  w <- N / sum(N)

  ## the idea is to return an efficiently stored symmetric matrix
  result <- list(genNrr=new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                            Dimnames=list(var.names, var.names),
                            x=rep(as.double(0), n.var*(n.var-1)/2+n.var)),
                 qOrders=qOrders)

  if (verbose && startTime["elapsed"] == 0)
    cat("Employing qOrders={", paste(paste(names(qOrders),
                                           qOrders, sep="="),
                                     collapse=", "),"}\n")

  elapsedTime <- 0
  for (idx in unique(datasetIdx)) {

    if (verbose && startTime["elapsed"] == 0)
      cat(sprintf("Data set %s\n", as.character(idx)))

    thisNrr <- qpgraph:::.qpNrr(X[datasetIdx == idx, ], qOrders[idx], I, restrict.Q,
                                fix.Q, nTests, alpha, pairup.i, pairup.j, verbose, identicalQs,
                                exact.test, use, tol, R.code.only, cl, startTime, nAdj2estimateTime)

    if (startTime["elapsed"] > 0) {
      elapsedTime <- elapsedTime + thisNrr["days"]*24*3600 + thisNrr["hours"]*3600 +
                     thisNrr["minutes"]*60 + thisNrr["seconds"]
      startTime <- proc.time()
    } else {
      ## this is necessary till we find out how to sum two dspMatrices getting a dspMatrix
      result[["genNrr"]] <- as(result[["genNrr"]] + w[idx] * thisNrr, "dspMatrix")

      if (return.all)
        result[[as.character(idx)]] <- thisNrr
    }
  }

  if (clusterSize > 1 && !is.null(cl)) {
    stopCl(cl)

    elapsedTime <- elapsedTime + (proc.time() - startTime)["elapsed"]
  }

  if (startTime["elapsed"] > 0) {
    d <- as.vector(floor(elapsedTime / (24*3600)))
    h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
    m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
    s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
    result <- c(days=d, hours=h, minutes=m, seconds=s)
  }

  return(result)
}



## function: qpEdgeNrr
## purpose: estimate the non-rejection rate for one pair of variables as the
##          fraction of tests that accept the null hypothesis of independence given
##          a set of randomly sampled q-order conditionals
## parameters: X - data set from where to estimate the non-rejection rate
##             i - index in S (row or column) matching one of the two variables
##             j - index in S (row or column) matching the other variable
##             q - partial-correlation order
##             I - indexes or names of the variables in X that are discrete
##             restrict.Q - indexes or names of variables to which the conditioning
##                          subsets Q should be restricted
##             fix.Q - indexes or names of variables that should be fixed within
##                     every conditioning subset Q
##             nTests - number of tests to perform
##             alpha - significance level of each test (Type-I error probability)
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             exact.test - employ an exact test when working with HMGMs
##             R.code.only - flag set to FALSE when using the C implementation
## return: an estimate of the non-rejection rate for the particular given pair of
##         variables

setGeneric("qpEdgeNrr", function(X, ...) standardGeneric("qpEdgeNrr"))

# X comes as an smlSet object
setMethod("qpEdgeNrr", signature(X="smlSet"),
          function(X, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                   nTests=100, alpha=0.05, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

            p <- as.integer(nrow(X))
            h <- as.integer(ncol(Biobase::pData(X)))
            sByChr <- sapply(GGBase::smList(X), ncol)
            cumsum_sByChr <- c(0, cumsum(sByChr))
            s <- sum(sByChr)
            n <- as.integer(ncol(X))
            fNames <- Biobase::featureNames(X)
            pNames <- colnames(Biobase::pData(X))
            sNamesByChr <- lapply(GGBase::smList(X), colnames)
            sNames <- unlist(sNamesByChr, use.names=FALSE)
            nLevels <- rep(NA, p+h)

            ## paste the involved phenotypic variables with the expression profiles
            ## just like when handling ExpressionSet objects and send this matrix
            ## with the smlSet object to the lower-level workhorse in order to
            ## minimize the memory footprint

            XP <- matrix(NA, nrow=ncol(X), ncol=0)
            I <-  NULL
            if (h > 0) { ## if there are phenotypic variables, they are allowed to
                         ## to be included in i, j or fix.Q
              if (is.character(i)) {
                if (!is.na(match(i, pNames))) { ## then 'i' refers to a phenotypic variable (cont. or discrete)
                  x <- Biobase::pData(X)[, i]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    x <- factor(x)
                    nLevels[p+ncol(XP)+1] <- nlevels(x)
                    I <- c(I, p+ncol(XP)+1)
                  }
                  XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, i)
                }
              } else if (i > p && i <= p+h) { ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i-p]
                cnames <- colnames(XP)
                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  x <- factor(x)
                  nLevels[p+ncol(XP)+1] <- nlevels(x)
                  I <- c(I, p+ncol(XP)+1)
                }
                XP <- cbind(XP, as.numeric(x))
                colnames(XP) <- c(cnames, pNames[i-p])
                i <- p+ncol(XP)
              } else if (i > p+h+s)
                stop(sprintf("i=%d is larger than the number of expression profiles (%d), phenotypic variables (%d) and genotypes together (%d+%d+%d=%d)\n", i, p, h, s, p+h+s))

              if (is.character(j)) {
                if (!is.na(match(j, pNames))) { ## then 'j' refers to a phenotypic variable (cont. or discrete)
                  x <- Biobase::pData(X)[, j]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    x <- factor(x)
                    nLevels[p+ncol(XP)+1] <- nlevels(x)
                    XP <- cbind(XP, as.numeric(x))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, j)
                }
              } else if (j > p && j <= p+h) { ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j-p]
                cnames <- colnames(XP)
                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  x <- factor(x)
                  nLevels[p+ncol(XP)+1] <- nlevels(x)
                  I <- c(I, p+ncol(XP)+1)
                }
                XP <- cbind(XP, as.numeric(x))
                colnames(XP) <- c(cnames, pNames[j-p])
                j <- p+ncol(XP)
              } else if (j > p+h+s)
                stop(sprintf("j=%d is larger than the number of expression profiles, phenotypic variables and genotypes together (%d+%d+%d=%d)\n", j, p, h, s, p+h+s))

              if (is.character(restrict.Q)) {
                mt <- match(restrict.Q, pNames)
                mt2 <- match(restrict.Q, colnames(XP)) ## avoid including a phen. var. twice
                for (k in mt[!is.na(mt) & is.na(mt2)]) {
                  x <- Biobase::pData(X)[, k]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    x <- factor(x)
                    nLevels[p+ncol(XP)+1] <- nlevels(x)
                    I <- c(I, p+ncol(XP)+1)
                  }
                  XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[k])
                }
              } else {
                for (k in which(restrict.Q > p & restrict.Q <= p+h)) {
                  if (is.na(match(pNames[restrict.Q[k]], colnames(XP)))) { ## avoid including a phen. var. twice
                    x <- Biobase::pData(X)[, restrict.Q[k]-p]
                    cnames <- colnames(XP)
                    if (is.character(x) || is.factor(x) || is.logical(x)) {
                      x <- factor(x)
                      nLevels[p+ncol(XP)+1] <- nlevels(x)
                      I <- c(I, p+ncol(XP)+1)
                    }
                    XP <- cbind(XP, as.numeric(x))
                    colnames(XP) <- c(cnames, pNames[restrict.Q[k]-p])
                    restrict.Q[k] <- p+ncol(XP)
                  }
                }
              }
              if (is.character(fix.Q)) {
                mt <- match(fix.Q, pNames)
                mt2 <- match(fix.Q, colnames(XP)) ## avoid including a phen. var. twice
                for (k in mt[!is.na(mt) & is.na(mt2)]) {
                  x <- Biobase::pData(X)[, k]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    x <- factor(x)
                    nLevels[p+ncol(XP)+1] <- nlevels(x)
                    I <- c(I, p+ncol(XP)+1)
                  }
                  XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[k])
                }
              } else {
                for (k in which(fix.Q > p & fix.Q <= p+h)) {
                  if (is.na(match(pNames[fix.Q[k]], colnames(XP)))) { ## avoid including a phen. var. twice
                    x <- Biobase::pData(X)[, fix.Q[k]-p]
                    cnames <- colnames(XP)
                    if (is.character(x) || is.factor(x) || is.logical(x)) {
                      x <- factor(x)
                      nLevels[p+ncol(XP)+1] <- nlevels(x)
                      I <- c(I, p+ncol(XP)+1)
                    }
                    XP <- cbind(XP, as.numeric(x))
                    colnames(XP) <- c(cnames, pNames[fix.Q[k]-p])
                    fix.Q[k] <- p+ncol(XP)
                  }
                }
              }
            } ## end if (h > 0)

            XEP <- t(Biobase::exprs(X))
            XEP <- cbind(XEP, XP)
            ph <- ncol(XEP) ## write down how many vars are expression profiles and phenotypes
            nLevels <- nLevels[1:ph]

            varNames <- c(colnames(XEP), sNames)
            Y <- 1:ph
            if (!is.null(I))
              Y <- Y[-I]

            param <- qpgraph:::.processParameters(varNames, ph, p+h, s, n, i=i, j=j, q=q,
                                                  I=I, Y=Y, restrict.Q=restrict.Q, fix.Q=fix.Q)
            i <- param$i
            j <- param$j
            I <- param$I
            Y <- param$Y
            restrict.Q <- param$restrict.Q
            fix.Q <- param$fix.Q

            if (length(intersect(restrict.Q, Y)) > 0)
              Y <- unique(c(intersect(i, Y), intersect(j, Y),
                            intersect(restrict.Q, Y), intersect(fix.Q, Y)))


            if (is.null(I) && !is.null(restrict.Q) &&
                all(c(i, j, restrict.Q, fix.Q) <= ph)) { ## only continuous variables are
                                                         ## involved in the calculations
              V <- c(i, j, setdiff(restrict.Q, c(i, j)), fix.Q)
              i <- 1L
              j <- 2L
              restrict.Q <- 2L+seq(along=setdiff(restrict.Q, c(i, j)))
              fix.Q <- 2L+length(restrict.Q)+seq(along=fix.Q)

              ## S <- qpCov(XEP[, V, drop=FALSE]) ## here is faster to calculate S for each margin
              S <- NULL

              qpgraph:::.qpEdgeNrr(XEP[, V, drop=FALSE], S, i, j, q, restrict.Q, fix.Q, nTests,
                                   alpha, R.code.only)
            } else {

              gLevels <- sum(unique(as.vector(as(GGBase::smList(X)[[1]][, 1:min(sByChr[1], 1000)], "matrix"))) > 0)

              qpgraph:::.qpEdgeNrrHMGMsml(GGBase::smList(X), cumsum_sByChr, s, gLevels,
                                          XEP, I, nLevels, Y, i, j, q, restrict.Q, fix.Q,
                                          nTests, alpha, exact.test, use, tol, R.code.only)
            }
          })

# X comes as an ExpressionSet object
setMethod("qpEdgeNrr", signature(X="ExpressionSet"),
          function(X, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                   nTests=100, alpha=0.05, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

            p <- as.integer(nrow(X))
            h <- as.integer(ncol(Biobase::pData(X)))
            n <- as.integer(ncol(X))
            fNames <- Biobase::featureNames(X)
            pNames <- colnames(Biobase::pData(X))

            XP <- matrix(NA, nrow=ncol(X), ncol=0)
            I <- NULL
            if (h > 0) { ## if there are phenotypic variables, they are allowed to
                         ## to be included in i, j or fix.Q
              if (is.character(i)) {
                if (!is.na(match(i, pNames))) { ## then 'i' refers to a phenotypic variable (cont. or discrete)
                  x <- Biobase::pData(X)[, i]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, i)
                }
              } else if (i > p && i <= p+h) { ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i-p]
                cnames <- colnames(XP)
                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                  I <- c(I, p+ncol(XP))
                } else
                  XP <- cbind(XP, as.numeric(x))
                colnames(XP) <- c(cnames, pNames[i-p])
                i <- p+ncol(XP)
              } else if (i > p+h)
                stop(sprintf("i=%d is larger than the number of expression profiles (%d) and phenotypic variables (%d) together (%d+%d=%d)\n", i, p, h, p+h))
              if (is.character(j)) {
                if (!is.na(match(j, pNames))) { ## then 'j' refers to a phenotypic variable (cont. or discrete)
                  x <- Biobase::pData(X)[, j]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, j)
                }
              } else if (j > p && j <= p+h) { ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j-p]
                cnames <- colnames(XP)
                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                  I <- c(I, p+ncol(XP))
                } else
                  XP <- cbind(XP, as.numeric(x))
                colnames(XP) <- c(cnames, pNames[j-p])
                j <- p+ncol(XP)
              } else if (j > p+h)
                stop(sprintf("j=%d is larger than the number of expression profiles (%d) and phenotypic variables (%d) together (%d+%d=%d)\n", j, p, h, p+h))
              if (is.character(restrict.Q)) {
                mt <- match(restrict.Q, pNames)
                mt2 <- match(restrict.Q, colnames(XP)) ## avoid including a phen. var. twice
                for (k in mt[!is.na(mt) & is.na(mt2)]) {
                  x <- Biobase::pData(X)[, k]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[k])
                }
              } else {
                for (k in which(restrict.Q > p)) {
                  if (is.na(match(pNames[restrict.Q[k]], colnames(XP)))) { ## avoid including a phen. var. twice
                    x <- Biobase::pData(X)[, restrict.Q[k]-p]
                    cnames <- colnames(XP)
                    if (is.character(x) || is.factor(x) || is.logical(x)) {
                      XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                      I <- c(I, p+ncol(XP))
                    } else
                      XP <- cbind(XP, as.numeric(x))
                    colnames(XP) <- c(cnames, pNames[restrict.Q[k]-p])
                    restrict.Q[k] <- p+ncol(XP)
                  }
                }
              }
              if (is.character(fix.Q)) {
                mt <- match(fix.Q, pNames)
                mt2 <- match(fix.Q, colnames(XP)) ## avoid including a phen. var. twice
                for (k in mt[!is.na(mt) & is.na(mt2)]) {
                  x <- Biobase::pData(X)[, k]
                  cnames <- colnames(XP)
                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                    I <- c(I, p+ncol(XP))
                  } else
                    XP <- cbind(XP, as.numeric(x))
                  colnames(XP) <- c(cnames, pNames[k])
                }
              } else {
                for (k in which(fix.Q > p)) {
                  if (is.na(match(pNames[fix.Q[k]], colnames(XP)))) { ## avoid including a phen. var. twice
                    x <- Biobase::pData(X)[, fix.Q[k]-p]
                    cnames <- colnames(XP)
                    if (is.character(x) || is.factor(x) || is.logical(x)) {
                      XP <- cbind(XP, as.numeric(factor(x, levels=unique(x))))
                      I <- c(I, p+ncol(XP))
                    } else
                      XP <- cbind(XP, as.numeric(x))
                    colnames(XP) <- c(cnames, pNames[fix.Q[k]-p])
                    fix.Q[k] <- p+ncol(XP)
                  }
                }
              }
            } ## end if (h > 0)

            X <- t(Biobase::exprs(X))
            X <- cbind(X, XP)
            varNames <- colnames(X)
            p <- ncol(X)

            if (is.null(I)) {
              param <- qpgraph:::.processParameters(varNames, p, p+h, 0, n, i=i, j=j, q=q,
                                                    restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              V <- 1:p
              if (!is.null(restrict.Q)) {
                V <- c(i, j, setdiff(restrict.Q, c(i, j)), fix.Q)
                i <- 1L
                j <- 2L
                restrict.Q <- 2L+seq(along=setdiff(restrict.Q, c(i, j)))
                fix.Q <- 2L+length(restrict.Q)+seq(along=fix.Q)
              }

              ## S <- qpCov(X[, V, drop=FALSE]) ## here is faster to calculate S for each margin
              S <- NULL

              qpgraph:::.qpEdgeNrr(X[, V, drop=FALSE], S, i, j, q, restrict.Q, fix.Q, nTests,
                                   alpha, R.code.only)
            } else {
              Y <- varNames
              if (is.character(I)) ## isn't it I at this point always integer ?
                Y <- setdiff(varNames, I)
              else
                Y <- (1:p)[-I]

              param <- qpgraph:::.processParameters(varNames, p, p+h, 0, n, i=i, j=j, q=q, I=I, Y=Y,
                                                      restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              I <- param$I
              Y <- param$Y
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              if (length(intersect(restrict.Q, Y)) > 0)
                Y <- unique(c(intersect(i, Y), intersect(j, Y),
                              intersect(restrict.Q, Y), intersect(fix.Q, Y)))

              nLevels <- rep(NA_integer_, times=p)
              nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

              ssd <- mapX2ssd <- NULL
              ## here is faster to calculate ssd for each margin
              ## missingMask <- apply(X, 2, function(x) any(is.na(x)))
              ## missingData <- any(missingMask)
              ## if (!missingData) {
              ##   ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
              ##   mapX2ssd <- match(varNames, colnames(ssd))
              ##   ## names(mapX2ssd) <- colnames(X) ## is this necessary ??
              ## }

              qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q,
                                       restrict.Q, fix.Q, nTests, alpha,
                                       exact.test, use, tol, R.code.only)
            }
          })

# X comes as a data frame
setMethod("qpEdgeNrr", signature(X="data.frame"),
          function(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                   nTests=100, alpha=0.05, long.dim.are.variables=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                   R.code.only=FALSE) {

            use <- match.arg(use)

            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)
            varNames <- colnames(X)
            p <- ncol(X)
            n <- nrow(X)

            if (is.null(I)) {
              param <- qpgraph:::.processParameters(varNames, p, p, 0, n, i=i, j=j, q=q,
                                                    restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              V <- 1:p
              if (!is.null(restrict.Q)) {
                V <- c(i, j, setdiff(restrict.Q, c(i, j)), fix.Q)
                i <- 1L
                j <- 2L
                restrict.Q <- 2L+seq(along=setdiff(restrict.Q, c(i, j)))
                fix.Q <- 2L+length(restrict.Q)+seq(along=fix.Q)
              }

              ## S <- qpCov(X[, V, drop=FALSE]) ## here is faster to calculate S for each margin
              S <- NULL

              qpgraph:::.qpEdgeNrr(X[, V, drop=FALSE], S, i, j, q, restrict.Q, fix.Q, nTests,
                                   alpha, R.code.only)
            } else {
              if (!is.character(I) && !is.numeric(I) && !is.integer(I))
                stop("I should be either variables names or indices\n")

              Y <- colnames(X)
              if (is.character(I))
                Y <- setdiff(colnames(X), I)
              else
                Y <- (1:ncol(X))[-I]

              param <- qpgraph:::.processParameters(varNames, p, p, 0, n, i=i, j=j, q=q, I=I, Y=Y,
                                                      restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              I <- param$I
              Y <- param$Y
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              if (length(intersect(restrict.Q, Y)) > 0)
                Y <- unique(c(intersect(i, Y), intersect(j, Y),
                              intersect(restrict.Q, Y), intersect(fix.Q, Y)))

              nLevels <- rep(NA_integer_, times=ncol(X))
              nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

              ssd <- mapX2ssd <- NULL
              ## here is faster to calculate ssd for each margin
              ## missingMask <- apply(X, 2, function(x) any(is.na(x)))
              ## missingData <- any(missingMask)
              ## if (!missingData) {
              ##   ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
              ##   mapX2ssd <- match(varNames, colnames(ssd))
              ##   ## names(mapX2ssd) <- colnames(X) ## is this necessary ??
              ## }

              qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q,
                                       restrict.Q, fix.Q, nTests, alpha,
                                       exact.test, use, tol, R.code.only)
            }
          })

          
# X comes as a matrix
setMethod("qpEdgeNrr", signature(X="matrix"),
          function(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                   nTests=100, alpha=0.05, long.dim.are.variables=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                   R.code.only=FALSE) {

            use <- match.arg(use)

            if (long.dim.are.variables &&
              sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)
            p <- ncol(X)
            varNames <- colnames(X)

            n <- nrow(X)
            if (is.null(I)) {
              param <- qpgraph:::.processParameters(varNames, p, p, 0, n, i=i, j=j, q=q,
                                                    restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              V <- 1:p
              if (!is.null(restrict.Q)) {
                V <- c(i, j, setdiff(restrict.Q, c(i, j)), fix.Q)
                i <- 1L
                j <- 2L
                restrict.Q <- 2L+seq(along=setdiff(restrict.Q, c(i, j)))
                fix.Q <- 2L+length(restrict.Q)+seq(along=fix.Q)
              }

              ## S <- qpCov(X[, V, drop=FALSE]) ## here is faster to calculate S for each margin
              S <- NULL

              qpgraph:::.qpEdgeNrr(X[, V, drop=FALSE], S, i, j, q, restrict.Q, fix.Q, nTests,
                                   alpha, R.code.only)
            } else {
              if (!is.character(I) && !is.numeric(I) && !is.integer(I))
                stop("I should be either variables names or indices\n")

              Y <- colnames(X)
              if (is.character(I))
                Y <- setdiff(colnames(X), I)
              else
                Y <- (1:ncol(X))[-I]

              param <- qpgraph:::.processParameters(varNames, p, p, 0, n, i=i, j=j, q=q, I=I, Y=Y,
                                                    restrict.Q=restrict.Q, fix.Q=fix.Q)
              i <- param$i
              j <- param$j
              I <- param$I
              Y <- param$Y
              restrict.Q <- param$restrict.Q
              fix.Q <- param$fix.Q

              if (length(intersect(restrict.Q, Y)) > 0)
                Y <- unique(c(intersect(i, Y), intersect(j, Y),
                              intersect(restrict.Q, Y), intersect(fix.Q, Y)))

              nLevels <- rep(NA_integer_, times=ncol(X))
              nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

              ssd <- mapX2ssd <- NULL
              ## here is faster to calculate ssd for each margin
              ## missingMask <- apply(X, 2, function(x) any(is.na(x)))
              ## missingData <- any(missingMask)
              ## if (!missingData) {
              ##   ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
              ##   mapX2ssd <- match(varNames, colnames(ssd))
              ##   ## names(mapX2ssd) <- colnames(X) ## is this necessary ??
              ## }

              qpgraph:::.qpEdgeNrrHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q,
                                       restrict.Q, fix.Q, nTests, alpha,
                                       exact.test, use, tol, R.code.only)
            }
          })

# X comes as an SsdMatrix (i.e., a sample covariance matrix calculated with qpCov())
setMethod("qpEdgeNrr", signature(X="SsdMatrix"),
          function(X, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                   nTests=100, alpha=0.05, R.code.only=FALSE) {

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)
            p <- ncol(X)
            varNames <- colnames(X)

            if (is.null(rownames(X)))
              rownames(X) <- colnames(X)

            param <- qpgraph:::.processParameters(varNames, p, p, 0, X@n, i=i, j=j, q=q,
                                                  restrict.Q=restrict.Q, fix.Q=fix.Q)
            i <- param$i
            j <- param$j
            restrict.Q <- param$restrict.Q
            fix.Q <- param$fix.Q

            qpgraph:::.qpEdgeNrr(NULL, X, i, j, q, restrict.Q, fix.Q, nTests,
                                 alpha, R.code.only)
          })

## ph contains the number of profile and phenotypic variables within varNames
## init_ph contains the initial number of profile and phenotypic variables on
## which integer values for i, j, restrict.Q or fix.Q could be based on.
.processParameters <- function(varNames, ph, init_ph, s, n, i, j, q, I=NULL, Y=NULL,
                               restrict.Q=NULL, fix.Q=NULL) {

  p <- length(varNames)

  if (q < 0)
    stop(paste("q=", q, " < 0"))

  if (q > n-3)
    stop(paste("q=", q, " > n-3=", n-3))

  if (is.character(i)) {
    if (is.na(match(i, varNames)))
      stop(sprintf("i=%s does not form part of the variable names of the data\n",i))
    i <- match(i, varNames)
  } else {
    if (i > init_ph && i <= init_ph+s)
      i <- i - init_ph + ph
  }

  if (is.character(j)) {
    if (is.na(match(j, varNames)))
      stop(sprintf("j=%s does not form part of the variable names of the data\n",j))
    j <- match(j, varNames)
  } else {
    if (j > init_ph && j <= init_ph+s)
      j <- j - init_ph + ph
  }

  if (!is.null(I)) {
    if (is.character(I)) {
      if (any(is.na(match(I, varNames))))
        stop("Some variables in I do not form part of the variable names of the data\n")
      I <- match(I, varNames)
    }
  }

  if (!is.null(Y)) {
    if (is.character(Y)) {
      if (any(is.na(match(Y, varNames))))
        stop("Some variables in Y do not form part of the variable names of the data\n")
      Y <- match(Y, varNames)
    }
  }

  if (!is.null(fix.Q)) {
    if (is.character(fix.Q)) {
      if (any(is.na(match(fix.Q, varNames))))
        stop("Some variables in fix.Q do not form part of the variable names of the data\n")
      fix.Q <- match(fix.Q, varNames)
    } else
      fix.Q[fix.Q > init_ph & fix.Q <= init_ph+s] <- fix.Q[fix.Q > init_ph & fix.Q <= init_ph+s] - init_ph + ph


    if (any(!is.na(match(c(i, j), fix.Q))))
      stop("The subset fix.Q cannot include any of the (i, j) variables.\n")

    if (q <= length(fix.Q))
      stop("q should be larger than the number of variables in fix.Q\n")
  }

  if (!is.null(restrict.Q)) {
    if (is.character(restrict.Q)) {
      if (any(is.na(match(restrict.Q, varNames))))
        stop("Some variables in restrict.Q do not form part of the variable names of the data\n")
      restrict.Q <- match(restrict.Q, varNames)
    } else
      restrict.Q[restrict.Q > init_ph & restrict.Q <= init_ph+s] <- restrict.Q[restrict.Q > init_ph & restrict.Q <= init_ph+s] - init_ph + ph
  }

  if (length(intersect(restrict.Q, fix.Q)) > 0)
    stop("The subsets restrict.Q and fix.Q should be disjoint.\n")

  return(list(i=i, j=j, I=I, Y=Y, restrict.Q=restrict.Q, fix.Q=fix.Q))
}

## IMPORTANT: .qpEdgeNrr() assumes that .processParameters() has been
##            previously called and all arguments related to variables
##            com as integers
.qpEdgeNrr <- function(X, S, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                       nTests=100, alpha=0.05, R.code.only=FALSE) {

  stopifnot(!is.null(X) || !is.null(S))

  if (!R.code.only) { ## assume restrict.Q and fix.Q are coordinately set!!!!
    return(qpgraph:::.qpFastEdgeNrr(X, S, i, j, q, restrict.Q, fix.Q, nTests, alpha));
  }

  Qm <- NA
  work.with.margin <- FALSE
  if (is.null(S)) {
    work.with.margin <- TRUE
    Qm <- seq(3, q+2)
  }

  p <- ifelse(!is.null(X), ncol(X), ncol(S))
  n <- ifelse(!is.null(X), nrow(X), S@n)
  V <- 1:p
  if (!is.null(restrict.Q))
    V <- restrict.Q

  if (!is.null(fix.Q))
    V <- setdiff(V, fix.Q)

  if (q > length(V)-2)
    stop(paste("q=", q, " > p-2=", p-2))

  q.fix <- length(fix.Q)
  mt <- match(c(i, j), V)
  mt <- mt[!is.na(mt)]
  if (length(mt) > 0)
    V <- V[-mt]

  thr <- qt(p=1-(alpha/2), df=n-q-2, lower.tail=TRUE, log.p=FALSE)
  lambda <- c()
  for (k in 1:nTests) {
    Q <- c(sample(V, q-q.fix, replace=FALSE), fix.Q)
    if (work.with.margin) {
      S <- qpgraph::qpCov(X[, c(i, j, Q)], corrected=TRUE)
      cit <- qpgraph:::.qpCItest(S, 1L, 2L, Qm, R.code.only=TRUE)
    } else
      cit <- qpgraph:::.qpCItest(S, as.integer(i), as.integer(j),
                                 as.integer(Q), R.code.only=TRUE)
    lambda  <- c(lambda, abs(cit$statistic))
  }

  nAcceptedTests <- sum(lambda < thr)

  return(nAcceptedTests / nTests)
}

## IMPORTANT: .qpEdgeNrrIdenticalQs() assumes that .processParameters() has been
##            previously called and all arguments related to variables
##            com as integers
## TODO: work with complete.obs on the margins
.qpEdgeNrrIdenticalQs <- function(S, Qs, S22invs, i=1, j=2, q=1, nTests=100, alpha=0.05,
                                  R.code.only=FALSE) {
  n <- S@n
  nActualTests <- 0
  thr    <- qt(p=1-(alpha/2), df=n-q-2, lower.tail=TRUE, log.p=FALSE)
  lambda <- c()
  for (k in 1:nTests) {
    if (sum(!is.na(match(c(i, j), Qs[[k]]))) == 0) { # those Q sets that include i or j are excluded
      Mmar    <- S[c(i, j, Qs[[k]]), c(i, j, Qs[[k]])]
      par.cov <- Mmar[1:2, 1:2] - Mmar[1:2, 3:(q+2)] %*% S22invs[[k]] %*% Mmar[3:(q+2), 1:2]
      par.cor <- cov2cor(par.cov)[1,2]
      t.value <- sqrt(n - q - 2) * par.cor / sqrt(1 - par.cor^2)
      lambda <- c(lambda, abs(t.value))
      nActualTests <- nActualTests + 1
    }
  }

  nAcceptedTests <- sum(lambda < thr)

  return(nAcceptedTests / nActualTests)
}

## IMPORTANT: .qpEdgeNrrHMGM() assumes that .processParameters() has been
##            previously called and all arguments related to variables
##            come as integers
.qpEdgeNrrHMGM <- function(X, I, nLevels, Y, ssdMat, mapX2ssdMat, i=1, j=2, q=1,
                           restrict.Q=NULL, fix.Q=NULL, nTests=100,
                           alpha=0.05, exact.test=TRUE, use=c("complete.obs", "em"),
                           tol=0.01, R.code.only=FALSE) {

  if (all(!is.na(match(c(i,j), I))))
    stop("i and j cannot be both discrete at the moment\n")

  if (!is.na(match(j, I))) { ## if any of (i,j) is discrete, it should be i
    tmp <- i
    i <- j
    j <- tmp
  }

  if (!R.code.only) {
    if (use == "em")
      stop("The EM algorithm is not yet implemented in the fast C version of the code for qpEdgeNrr(). Please set either use=\"complete.obs\" (default value) or R.code.only=TRUE\n")

    return(qpgraph:::.qpFastEdgeNrrHMGM(X, I, nLevels, Y, ssdMat, mapX2ssdMat, i, j, q,
                                        restrict.Q, fix.Q, nTests, alpha, exact.test))
  }

  p <- ncol(X)
  V <- 1:p
  if (!is.null(restrict.Q))
    V <- restrict.Q

  if (!is.null(fix.Q))
    V <- setdiff(V, fix.Q)

  if (q > length(V)-2)
    stop(paste("q=", q, " > p-2=", p-2))

  q.fix <- length(fix.Q)
  mt <- match(c(i, j), V)
  mt <- mt[!is.na(mt)]
  if (length(mt) > 0)
    V <- V[-mt]

  problematicQ <- NA
  nActualTests <- 0
  lambda <- a <- b <- thr <- rep(NA, times=nTests)
  for (k in 1:nTests) {
    Q <- c(sample(V, q-q.fix, replace=FALSE), fix.Q)
    cit <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssdMat, mapX2ssdMat, as.integer(i),
                                   as.integer(j), as.integer(Q), exact.test,
                                   use, tol, R.code.only=TRUE)
    if (!is.nan(cit$statistic)) {
      lambda[k] <- cit$statistic
      if (exact.test) {
        a[k] <- cit$parameter["a"]
        b[k] <- cit$parameter["b"]
        if (k > 1 && a[k] == a[k-1] && b[k] == b[k-1])
          thr[k] <- thr[k-1]
        else
          thr[k] <- qbeta(p=alpha, shape1=a[k], shape2=b[k], lower.tail=TRUE)
      }
      nActualTests <- nActualTests + 1
    } else
      problematicQ <- Q
  }

  nAcceptedTests <- NA
  if (exact.test) {
    nAcceptedTests <- sum(lambda > thr, na.rm=TRUE)
  } else {
    thr <- qchisq(p=(1-alpha), df=1, lower.tail=TRUE)
    nAcceptedTests <- sum(lambda < thr, na.rm=TRUE)
  }


  if (nActualTests < nTests)
    warning(paste(sprintf("Non-rejection rate estimation between i=%s and j=%s with q=%d was based on %d out of %d requested tests.\n",
                          colnames(X)[i], colnames(X)[j], q, nActualTests, nTests),
                  sprintf("For instance, the CI test between i=%s and j=%s given Q={",
                          colnames(X)[i], colnames(X)[j]),
                  paste(colnames(X)[problematicQ], collapse=", "),
                  "}, could not be calculated. Try with smaller Q subsets or increase n if you can.\n",
                  sep=""))

  return(nAcceptedTests / nActualTests)
}

## automatic calculation of genotype levels
## sum(unique(as.vector(as(smList(c20)[[1]][, 1:1000], "matrix"))) > 0)

## IMPORTANT: .qpEdgeNrrHMGMsml() assumes that .processParameters() has been
##            previously called and all arguments related to variables
##            com as integers
.qpEdgeNrrHMGMsml <- function(X, cumsum_sByChr, s, gLevels, XEP, I, nLevels, Y,
                              i=1L, j=2L, q=1L, restrict.Q=NULL, fix.Q=NULL,
                              nTests=100, alpha=0.05, exact.test=TRUE,
                              use=c("complete.obs", "em"), tol=0.01,
                              R.code.only=FALSE) {
  n <- nrow(XEP)
  ph <- ncol(XEP)
  p <- ph + s

  iDiscrete <- !is.na(match(i, I)) | i > ph
  jDiscrete <- !is.na(match(j, I)) | j > ph
  if (iDiscrete && jDiscrete)
    stop("i and j cannot be both discrete at the moment\n")

  ssd <- mapY2ssd <- NULL
  ## ssd <- qpCov(XEP[, Y, drop=FALSE], corrected=FALSE)
  ## ## we need to map smlSet-level Y indices to the ssd dimension
  ## mapY2ssd <- rep(NA, length=ph)
  ## mapY2ssd[Y] <- 1:ncol(ssd)

  if (!R.code.only) {
    return(qpgraph:::.qpFastEdgeNrrHMGMsml(X, cumsum_sByChr, s, gLevels, XEP, I, nLevels,
                                           Y, ssd, mapY2ssd, i, j, q, restrict.Q,
                                           fix.Q, nTests, alpha, exact.test))
  }

  if (!is.na(match(j, I)) || j > ph) { ## if any of (i,j) is discrete, it should be i
    tmp <- i
    i <- j
    j <- tmp
    ## tmp <- nLevels[i]
    ## nLevels[i] <- nLevels[j]
    ## nLevels[j] <- tmp
  }

  V <- 1:p
  if (!is.null(restrict.Q))
    V <- restrict.Q

  if (!is.null(fix.Q))
    V <- setdiff(V, fix.Q)

  if (q > length(V)-2)
    stop(paste("q=", q, " > p-2=", p-2))

  q.fix <- length(fix.Q)
  mt <- match(c(i, j), V)
  mt <- mt[!is.na(mt)]
  if (length(mt) > 0)
    V <- V[-mt]

  Xsub <- matrix(NA_real_, nrow=n, ncol=q+2)
  nLevelsSub <- rep(NA, q+2)
  Yij <- Iij <- NULL
  if (i > ph) {
    selChr <- sum(cumsum_sByChr < i-ph)
    x <- as(X[[selChr]][, i-ph-cumsum_sByChr[selChr]], "numeric")[ ,1]+1
    if (any(x > 3))
      warning(sprintf("i=%s has uncertain genotype calls which are treated here as missing", i))
    x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
    Xsub[, 1] <- x
    Iij <- 1
    nLevelsSub[1] <- gLevels
  } else {
    Xsub[, 1] <- XEP[, i]
    Yij <- i
  }
  i <- 1L 

  if (j > ph) {
    selChr <- sum(cumsum_sByChr < j-ph)
    x <- as(X[[selChr]][, j-ph-cumsum_sByChr[selChr]], "numeric")[ ,1]+1
    if (any(x > 3))
      warning(sprintf("j=%s has uncertain genotype calls which are treated here as missing", j))
    x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
    Xsub[, 2] <- x
    Iij <- 1 ## assuming only i or j can be discrete, not both at the same time
    nLevelsSub[2] <- gLevels
  } else {
    Xsub[, 2] <- XEP[, j]
    Yij <- j
  }
  j <- 2L

  problematicQ <- NA
  nActualTests <- 0
  lambda <- a <- b <- thr <- rep(NA, times=nTests)
  for (k in 1:nTests) {
    Q <- c(sample(V, q-q.fix, replace=FALSE), fix.Q)
    w <- which(Q <= ph)
    nijep <- length(w)
    Xsub[, 2L+seq(along=w)] <- XEP[, Q[w]]
    Isub <- c(Iij, 2L+seq(along=w)[na.omit(match(I, Q))])
    nLevelsSub[2L+seq(along=w)[na.omit(match(I, Q))]] <- nLevels[Q[na.omit(match(I, Q))]]
    Qw <- Q[which(Q > ph)]
    for (m in seq(along=Qw)) {
      selChr <- sum(cumsum_sByChr < Qw[m]-ph)
      x <- as(X[[selChr]][, Qw[m]-ph-cumsum_sByChr[selChr]], "numeric")[, 1]+1
      x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
      Xsub[, 2+nijep+m] <- x
      Isub <- c(Isub, 2+nijep+m)
      nLevelsSub[2+nijep+m] <- gLevels
    }
 
    ## map2ssd maps columns of continuous variables in Xsub to columns in ssd
    ## map2ssd <- rep(NA, length=q+2)
    Ysub <- (1:(q+2))
    if (length(Isub) > 0)
      Ysub <- Ysub[-Isub]
    ## map2ssd[Ysub] <- mapY2ssd[intersect(c(Yij, Q), Y)] ## map Xsub coord. to ssd coord.

    Q <- 2L+seq(along=Q)

    cit <- qpgraph:::.qpCItestHMGM(Xsub, Isub, nLevelsSub, Ysub, ssd, mapY2ssd,
                                   i, j, Q, exact.test, use, tol, R.code.only=TRUE)
    if (!is.nan(cit$statistic)) {
      lambda[k] <- cit$statistic
      if (exact.test) {
        a[k] <- cit$parameter["a"]
        b[k] <- cit$parameter["b"]
        if (k > 1 && a[k] == a[k-1] && b[k] == b[k-1])
          thr[k] <- thr[k-1]
        else
          thr[k] <- qbeta(p=alpha, shape1=a[k], shape2=b[k], lower.tail=TRUE)
      }
      nActualTests <- nActualTests + 1
    } else
      problematicQ <- Q
  }

  nAcceptedTests <- NA
  if (exact.test) {
    nAcceptedTests <- sum(lambda > thr, na.rm=TRUE)
  } else {
    thr <- qchisq(p=(1-alpha), df=1, lower.tail=TRUE)
    nAcceptedTests <- sum(lambda < thr, na.rm=TRUE)
  }


  if (nActualTests < nTests)
    ## UPDATE HOW VARIABLE NAMES ARE HERE RETRIEVED !!!
    warning(paste(sprintf("Non-rejection rate estimation between i=%s and j=%s with q=%d was based on %d out of %d requested tests.\n",
                          colnames(X)[i], colnames(X)[j], q, nActualTests, nTests),
                  sprintf("For instance, the CI test between i=%s and j=%s given Q={",
                          colnames(X)[i], colnames(X)[j]),
                  paste(colnames(X)[problematicQ], collapse=", "),
                  "}, could not be calculated. Try with smaller Q subsets or increase n if you can.\n",
                  sep=""))

  return(nAcceptedTests / nActualTests)
}



## function: qpHist
## purpose: plot the distribution of non-rejection rates
## parameters: nrrMatrix - matrix of non-rejection rates
##             A - adjacency matrix whose present and missing edges will be employed
##                 to show separately the distribution of non-rejection rates
##             freq - logical; if TRUE, the histograms show frequencies (counts)
##                    of occurrence of the different non-rejection rate values;
##                    if FALSE, then probability densities are plotted
## return: none

qpHist <- function(nrrMatrix, A=NULL,
                   titlehist = "all estimated\nnon-rejection rates", freq=TRUE) {
  # all
  nrr <- nrrMatrix[upper.tri(nrrMatrix)]
  nrr_rg <- range(nrr)
  if(is.null(A)){
    hist(nrr, col="yellow", main=titlehist, xlab="non-rejection rate", freq=freq)
  } else {
    # only beta
    T <- as.matrix(nrrMatrix) ## till [<- works with dspMatrix matrices
    T[!as.matrix(A)] <- -1
    xbeta <- T[upper.tri(T)]
    xbeta <- xbeta[xbeta != -1]
    # not beta
    T <- as.matrix(nrrMatrix) ## till [<- works with dspMatrix matrices
    T[as.matrix(A)] <- -1
    xnotbeta <- T[upper.tri(T)]
    xnotbeta <- xnotbeta[xnotbeta != -1]
    # plots
    split.screen(c(2, 2))
    screen(1)
    H <- hist(nrr, plot=FALSE)
    if(freq){
       yl <- c(0, max(H$counts))
    }else{
       yl < NULL
    }
    hist(nrr, freq=freq, col="yellow", xlim=nrr_rg, ylim=yl, main=titlehist,
         xlab="non-rejection rates")
    screen(2)
    boxplot(xbeta, xnotbeta, col=c("red", "skyblue"), ylab="non-rejection rate",
            naxt="n")
    axis(1,at=1,"present\nedges")
    axis(1,at=2,"missing\nedges")
    screen(3)
    hist(xbeta, freq=freq,  col="yellow", xlim=nrr_rg, ylim=yl,
         main="present edges\nnon-rejection rates", xlab="non-rejection rates",
         breaks=length(H$breaks))
    screen(4)
    hist(xnotbeta, freq=freq, col="yellow", xlim=nrr_rg, ylim=yl,
         main="missing edges\nnon-rejection rates", xlab="non-rejection rates",
         breaks=length(H$breaks))
    close.screen(all.screens=TRUE)
  }
}



## function: qpGraph
## purpose: obtain a qp-graph from a matrix of non-rejection rates
## parameters: nrrMatrix - matrix of non-rejection rates
##             threshold - threshold for edge removal
##             topPairs - form the qp-graph with a number topPairs of edges
##                        starting from the top of the ranking defined by
##                        the values in nrrMatrix
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             return.type - type of data structure on which the graph
##                           should be returned, either an adjacency matrix,
##                           a matrix with the list of edges, a graphNEL object
##                           a graphAM object or a graphBAM object
## return: adjacency matrix of the qp-graph

qpGraph <- function(nrrMatrix, threshold=NULL, topPairs=NULL, pairup.i=NULL,
                    pairup.j=NULL, return.type=c("adjacency.matrix", "edge.list",
                    "graphNEL", "graphAM", "graphBAM")) {

  return.type <- match.arg(return.type)

  ## by now we need to coerce the dspMatrix into a regular matrix
  ## hopefully this can be removed in the near future
  nrrMatrix <- as(nrrMatrix, "matrix")

  n.var <- nrow(nrrMatrix)

  if (is.null(colnames(nrrMatrix))) {
    vertex.labels <- as.character(1:n.var)
  } else {
    vertex.labels <- colnames(nrrMatrix)
  }

  if (is.null(threshold) && is.null(topPairs))
    stop("either threshold or topPairs should be set different to NULL\n")

  if (!is.null(threshold) && !is.null(topPairs))
    stop("only either threshold or topPairs can be set different to NULL\n")

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(nrrMatrix)))
      stop("when using pairup.i and pairup.j, nrrMatrix should have row and column names\n")

    var.names <- colnames(nrrMatrix)
    pairup.i <- match(pairup.i, var.names)
    if (sum(is.na(pairup.i)) > 0)
      stop("pairup.i is not a subset of the variables forming the data\n")
    pairup.j <- match(pairup.j, var.names)
    if (sum(is.na(pairup.j)) > 0)
      stop("pairup.j is not a subset of the variables forming the data\n")

    pairup.ij.int <- intersect(pairup.i, pairup.j)
    pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
    pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

    nomeasurementsMask <- matrix(FALSE,nrow=n.var,ncol=n.var)
    nomeasurementsMask[as.matrix(
                       expand.grid(pairup.ij.int,
                                   c(pairup.i.noint, pairup.j.noint)))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.i.noint, pairup.j.noint))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.ij.int, pairup.ij.int))] <- TRUE
    diag(nomeasurementsMask) <- FALSE
    nomeasurementsMask <- nomeasurementsMask | t(nomeasurementsMask)
    nomeasurementsMask <- !nomeasurementsMask
    nrrMatrix[nomeasurementsMask] <- NA
  }

  # not-available NRRs imply no edges
  nrrMatrix[is.na(nrrMatrix)] <- 1.0

  if (!is.null(threshold)) {
    A <- nrrMatrix <= threshold
  } else { # topPairs
    nrrUppTriMatrix <- nrrMatrix[upper.tri(nrrMatrix)]
    rowUppTri <- row(nrrMatrix)[upper.tri(nrrMatrix)]
    colUppTri <- col(nrrMatrix)[upper.tri(nrrMatrix)]
    orderedMeasurementsIdx <- sort(nrrUppTriMatrix, index.return=TRUE,
                                   decreasing=FALSE)$ix
    ranking <- cbind(rowUppTri[orderedMeasurementsIdx],
                     colUppTri[orderedMeasurementsIdx])
    A <- matrix(FALSE, nrow=n.var, ncol=n.var)
    A[ranking[1:topPairs, ]] <- TRUE
    A[cbind(ranking[1:topPairs, 2], ranking[1:topPairs, 1])] <- TRUE
  }

  rownames(A) <- colnames(A) <- vertex.labels
  diag(A) <- FALSE # whatever the threshold is the graph should have no loops

  if (return.type == "adjacency.matrix") {
    return(Matrix::Matrix(A))
  } else if (return.type == "edge.list") {
    m <- cbind(vertex.labels[row(A)[upper.tri(A) & A]], vertex.labels[col(A)[upper.tri(A) & A]])
    colnames(m) <- c("i", "j")
    return(m)
  } else {
    from <- vertex.labels[row(A)[upper.tri(A) & A]]
    to <- vertex.labels[col(A)[upper.tri(A) & A]]
    df <- data.frame(from=from, to=to, weight=rep(1, length(from)), stringsAsFactors=FALSE)
    g <- graph::graphBAM(df, nodes=vertex.labels)
    if (return.type == "graphNEL") {
      ## m <- cbind(vertex.labels[row(A)[upper.tri(A) & A]], vertex.labels[col(A)[upper.tri(A) & A]])
      ## m <- rbind(m, m[, 2:1])
      ## edL <- lapply(split(m[,1], m[,2]), unique)
      ## g <- new("graphNEL", nodes=names(edL), edgeL=edL, edgemode="undirected")
      g <- as(g, "graphNEL")
      return(g)
    } else if (return.type == "graphAM") {
      ## g <- new("graphAM", adjMat=A+0, edgemode="undirected", values=list(weight=1))
      g <- as(g, "graphAM")
      return(g)
    }
    return(g) ## graphBAM
  }

  return(NA)
}



## function: qpAnyGraph
## purpose: obtain a graph from a matrix of measurements by thresholding on
##          these measurements
## parameters: measurementsMatrix - matrix of pairwise associations
##             threshold - threshold for edge removal
##             remove - direction of the removal with the threshold
##             topPairs - form the graph with a number topPairs of edges
##                        starting from the top of the ranking defined by
##                        the measurementsMatrix
##             decreasing - logical, only applies when topPairs is set; if TRUE
##                          then the ranking is made in decreasing order; if
##                          FALSE then is made in increasing order
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             return.type - type of data structure on which the graph
##                           should be returned, either an adjacency matrix,
##                           a matrix with the list of edges, a graphNEL object,
##                           a graphAM object or a graphBAM object
## return: adjacency matrix of the qp-graph

qpAnyGraph <- function(measurementsMatrix, threshold=NULL, remove=c("below", "above"),
                       topPairs=NULL, decreasing=TRUE, pairup.i=NULL, pairup.j=NULL,
                       return.type=c("adjacency.matrix", "edge.list", "graphNEL",
                                     "graphAM", "graphBAM")) {

  remove <- match.arg(remove)
  return.type <- match.arg(return.type)

  ## by now we need to coerce the dspMatrix into a regular matrix
  ## hopefully in the near future we can do also [<- on dspMatrix matrices
  measurementsMatrix <- as(measurementsMatrix, "matrix")

  n.var <- nrow(measurementsMatrix)

  if (is.null(colnames(measurementsMatrix))) {
    vertex.labels <- as.character(1:n.var)
  } else {
    vertex.labels <- colnames(measurementsMatrix)
  }

  if (is.null(threshold) && is.null(topPairs))
    stop("either threshold or topPairs should be set different to NULL\n")

  if (!is.null(threshold) && !is.null(topPairs))
    stop("only either threshold or topPairs can be set different to NULL\n")

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(measurementsMatrix)))
      stop("when using pairup.i and pairup.j, measurementsMatrix should have row and column names\n")

    var.names <- colnames(measurementsMatrix)
    pairup.i <- match(pairup.i, var.names)
    if (sum(is.na(pairup.i)) > 0)
      stop("pairup.i is not a subset of the variables forming the data\n")
    pairup.j <- match(pairup.j, var.names)
    if (sum(is.na(pairup.j)) > 0)
      stop("pairup.j is not a subset of the variables forming the data\n")

    pairup.ij.int <- intersect(pairup.i, pairup.j)
    pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
    pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

    nomeasurementsMask <- matrix(FALSE,nrow=n.var,ncol=n.var)
    nomeasurementsMask[as.matrix(
                       expand.grid(pairup.ij.int,
                                   c(pairup.i.noint, pairup.j.noint)))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.i.noint, pairup.j.noint))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.ij.int, pairup.ij.int))] <- TRUE
    diag(nomeasurementsMask) <- FALSE
    nomeasurementsMask <- nomeasurementsMask | t(nomeasurementsMask)
    nomeasurementsMask <- !nomeasurementsMask
    measurementsMatrix[nomeasurementsMask] <- NA
  }

  ## non-available measurements imply no edges
  measurementsMatrix[is.na(measurementsMatrix)] <- NA

  if (!is.null(threshold)) {
    if (remove == "below")
      A <- measurementsMatrix >= threshold
    else
      A <- measurementsMatrix <= threshold
  } else { ## topPairs
    if (decreasing)
      bottomValue <- min(measurementsMatrix,na.rm=TRUE) - 1
    else
      bottomValue <- max(measurementsMatrix,na.rm=TRUE) + 1

    measurementsMatrix[is.na(measurementsMatrix)] <- bottomValue
    measurementsUppTriMatrix <- measurementsMatrix[upper.tri(measurementsMatrix)]
    rowUppTri <- row(measurementsMatrix)[upper.tri(measurementsMatrix)]
    colUppTri <- col(measurementsMatrix)[upper.tri(measurementsMatrix)]
    orderedMeasurementsIdx <- sort(measurementsUppTriMatrix, index.return=TRUE,
                                   decreasing=decreasing)$ix
    ranking <- cbind(rowUppTri[orderedMeasurementsIdx],
                     colUppTri[orderedMeasurementsIdx])
    A <- matrix(FALSE, nrow=n.var, ncol=n.var)
    A[ranking[1:topPairs,]] <- TRUE
    A[cbind(ranking[1:topPairs,2], ranking[1:topPairs,1])] <- TRUE
  }

  A[is.na(A)] <- FALSE
  rownames(A) <- colnames(A) <- vertex.labels
  diag(A) <- FALSE ## whatever the threshold is the graph should have no loops

  if (return.type == "adjacency.matrix") {
    ## require(Matrix)
    return(Matrix(A))
  } else if (return.type == "edge.list") {
    m <- cbind(row(A)[upper.tri(A) & A],col(A)[upper.tri(A) & A])
    colnames(m) <- c("i","j")
    return(m)
  } else if (return.type == "graphNEL") {
    vertices <- unique(c(vertex.labels[row(A)[upper.tri(A) & A]], vertex.labels[col(A)[upper.tri(A) & A]]))
    edL <- vector("list", length=length(vertices))
    names(edL) <- vertices
    for (v in vertices)
      edL[[v]] <- list(edges=vertex.labels[A[v, ]], weights=rep(1, sum(A[v, ])))
    g <- new("graphNEL",nodes=vertices, edgeL=edL, edgemode="undirected")
    return(g)
  } else if (return.type == "graphAM") {
    ## require(graph)
    g <- new("graphAM", adjMat=A+0, edgemode="undirected", values=list(weight=1))
    return(g)
  } else if (return.type == "graphAM") {
    vertices <- unique(c(vertex.labels[row(A)[upper.tri(A) & A]], vertex.labels[col(A)[upper.tri(A) & A]]))
    edL <- vector("list", length=length(vertices))
    names(edL) <- vertices
    for (v in vertices)
      edL[[v]] <- list(edges=vertex.labels[A[v, ]], weights=rep(1, sum(A[v, ])))
    g <- new("graphNEL",nodes=vertices, edgeL=edL, edgemode="undirected")
    return(as(g, "graphBAM"))
  }

  return(NA)
}



## function: qpGraphDensity
## purpose: calculate and plot the graph density as function of the non-rejection
##          rate
## parameters: nrrMatrix - matrix of non-rejection rates
##             threshold.lim - range of threshold values
##             breaks - either a number of threshold bins or a vector of
##             threshold breakpoints
##             plot - flag setting to make a plot of the result
##             qpGraphDensityOutput - output from a previous call to
##                                    qpGraphDensity, this allows one to
##                                    plot the result changing some of the
##                                    plotting parameters without having to do
##                                    the calculation again
##             density.digits - number of digits in the reported graph densities
##             titlegd - title to be shown in the plot
## return: a list with the graph density as function of threshold and an estimate
##         of the sparseness of the resulting qp-graphs across the thresholds

qpGraphDensity <- function(nrrMatrix, threshold.lim=c(0,1), breaks=5,
                           plot=TRUE, qpGraphDensityOutput=NULL,
                           density.digits=0,
                           titlegd="graph density as function of threshold") {

  if (is.null(qpGraphDensityOutput)) {

    if (length(breaks) == 1) {
      len <- threshold.lim[2] - threshold.lim[1]
      br <- seq(threshold.lim[1],threshold.lim[2],by=len/breaks)
    } else {
      br <- breaks
      threshold.lim = range(br)
    }

    matgdthr <- matrix(rep(0,length(br)*2),nrow=length(br),ncol=2)
    colnames(matgdthr) <- c("density", "threshold")
    n.var <- nrow(nrrMatrix)
    n.adj <- n.var*(n.var-1)/2

    ## by now we coerce this to a regular matrix
    nrrMatrix <- as(nrrMatrix, "matrix")

    for (i in 1:length(br)) {
      threshold <- br[i]
      nrrMatrix[is.na(nrrMatrix)] <- 1.0 # non-available NRRs imply no edges
      A <- nrrMatrix <= threshold
      diag(A) <- FALSE # if the threshold is 1.0 the resulting qp-graph
                       # will be the complete undirected graph but
                       # still it should have no loops
      n.edg <- sum(A) / 2
      gd <- (n.edg / n.adj) * 100
      matgdthr[i,] <- c(gd,threshold)
    }

  } else {
    br <- qpGraphDensityOutput$data[,2]
    matgdthr <- qpGraphDensityOutput$data
  }

  linetype <- 1
  label <- "graph density"

  if (plot == TRUE) {
    plot(br,matgdthr[,1],type="o",xlim=threshold.lim,
         ylim=range(1,matgdthr[,1]),lwd=2,lty=linetype,
         xlab="threshold",ylab="graph density",main=titlegd)
    pct <- round(matgdthr[,1],digits=density.digits)
    text(br,matgdthr[,1],lab=paste(pct,"%",sep=""),pos=1,cex=.7)
    legend(min(threshold.lim),max(matgdthr[,1]),label,lty=linetype,lwd=2,pch=1)
  }

  m <- matgdthr[,c(2,1)]
  m[,2] <- m[,2] / 100
  f <- approxfun(m)
  area <- integrate(f,min(m[,1]),max(m[,1]))

  invisible(list(data=matgdthr,sparseness=1-area$value))
}



## function: qpCliqueNumber
## purpose: calculate the size of the largest maximal clique in a given undirected graph
## parameters: g - either a graphNEL, graphAM, graphBAM object or an adjacency matrix of the graph
##             exact.calculation - flag that when set to TRUE the exact maximum clique
##                                 size is calculated and when set to FALSE a lower
##                                 bound is calculated instead
##             return.vertices - returns one set of vertices forming a maximal clique of the
##                               maximum size when this flag is set to TRUE
##             approx.iter - number of iterations performed to calculate
##                           the lower bound on the clique number of
##                           each graph (exact.calculation is FALSE)
##             verbose - show progress on the clique number calculation
##             R.code.only - flag set to FALSE when using the C implementation
## return: the size of the largest maximal clique in the given graph, also known as
##         its clique number

qpCliqueNumber <- function(g, exact.calculation=TRUE, return.vertices=FALSE,
                           approx.iter=100, verbose=TRUE, R.code.only=FALSE) {

  if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    ## require(graph)
    if (graph::edgemode(g) != "undirected")
      stop("g should be an undirected graph\n")

    A <- as(g, "matrix") == 1
  } else if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    A <- g
    p <- (d <- dim(A))[1]
    if (p != d[2])
      stop("g is not an squared adjacency matrix nor a graphNEL, graphAM or graphBAM object\n")

    if (!isSymmetric(A))
      stop("g is not a symmetric adjacency matrix nor a graphNEL, graphAM or graphBAM object\n")
  } else
    stop("g should be either a graphNEL object, a graphAM object, a graphBAM object or a boolean adjacency matrix\n")

  if (exact.calculation && R.code.only)
    stop("R code is only available for the lower bound approximation and not for the exact calculation\n");

  n.var <- nrow(A)
  n.possibleedges <- (n.var * (n.var-1)) / 2

  if (!any(A)) {
    maximum_clique <- 1
    if (return.vertices) {
        maximum_clique <- list(size=clique.number,vertices=1)
    }
    return(maximum_clique)
  }

  if (sum(A[upper.tri(A)]) == n.possibleedges) {
    maximum_clique <- n.var
    if (return.vertices) {
        maximum_clique <- list(size=clique.number,vertices=1:n.var)
    }
    return(maximum_clique)
  }

  maximum_clique <- 0

  if (exact.calculation == TRUE) {

    A <- A == 1 ## make sure we get a boolean matrix

    maximum_clique <-
      qpgraph:::.qpCliqueNumberOstergard(as.matrix(A),
                                         return.vertices=return.vertices,
                                         verbose=verbose)
  } else {
    if (!R.code.only)
      maximum_clique <-
        qpgraph:::.qpCliqueNumberLowerBound(as.matrix(A),
                                            return.vertices=return.vertices,
                                            approx.iter=approx.iter,
                                            verbose=verbose)
    else {
      if (verbose) {
        cat("calculating lower bound on the maximum clique size\n")
      }

      clique.number <- 0
      clique.vertices <- c()

      A <- as.matrix(A) + 0 ## make sure we get a 0-1 matrix
      deg <- sort(rowsum(A, rep(1,n.var)), index.return=TRUE,
                  decreasing=TRUE) ## order by degree

      ppct <- -1
      for (i in 1:approx.iter) {

        pdeg <- deg$ix
        if (i %% n.var + 1 > 1) {
          sset <- sample(1:n.var, i %% n.var + 1, replace=FALSE) ## we alter the order of the ranking
          ssetelem <- pdeg[sset]                                 ## by degree with increasing levels
          ssetelem <- sample(ssetelem)                           ## of randomness cyclically
          pdeg[sset] <- ssetelem
        }
        clq <- c(pdeg[1])
        j <- 2
        for (j in 2:n.var) {
          v <- pdeg[j]
          clq2 <- c(clq,v)
          if (sum(A[clq2,clq2]) == length(clq2)*length(clq2)-length(clq2)) {
            clq <- clq2
          }
        }
                                                                                            
        if (length(clq) > clique.number) {
          clique.number <- length(clq)
          clique.vertices <- clq
        }

        if (verbose) {
          pct <- floor((i*100)/approx.iter)
          if (pct != ppct) {
            if (pct %% 10 == 0) {
              cat(pct)
            } else {
              cat(".")
            }
            ppct = pct
          }
        }
      }

      if (verbose) {
        cat("\n")
      }

      maximum_clique <- clique.number
      if (return.vertices) {
        maximum_clique <- list(size=clique.number,vertices=clique.vertices)
      }
    }

  }

  return(maximum_clique)
}



## function: qpClique
## purpose: calculate and plot the maximum clique size as function of the
##          non-rejection rate
## parameters: nrrMatrix - matrix of non-rejection rates
##             n - sample size
##             threshold.lim - range of threshold values
##             breaks - either a number of threshold bins or a vector of
##                      threshold breakpoints
##             plot - flag that when set it makes a plot of the result
##             exact.calculation - flag that when set to TRUE the exact maximum
##                                 clique size is calculated and when set to
##                                 FALSE a lower bound is calculated instead
##             approx.iter - number of iterations performed to
##                           calculate the lower bound on the
##                           maximum clique size of each graph
##                           (i.e., only applies when
##                           exact.calculation is FALSE)
##             qpCliqueOutput - output from a previous call to qpClique, this
##                              allows one to plot the result changing some of
##                              the plotting parameters without having to do
##                              the calculation again
##             density.digits - number of digits in the reported graph densities
##             logscale.clqsize - whether the maximum clique size in the y-axis
##                                should be plotted in logarithmic scale
##             titleclq - title to be shown in the plot
##             verbose - show progress when doing the calculation
## return: a list with the maximum clique size and the graph density as function
##         of the threshold, an estimate of the complexity of the resulting
##         qp-graphs across the thresholds, the threshold on the non-rejection
##         rate that provides a maximum clique size strictly smaller than the
##         sample size n, and the resulting maximum clique size

qpClique <- function(nrrMatrix, n=NA, threshold.lim=c(0,1), breaks=5, plot=TRUE,
                     exact.calculation=TRUE, approx.iter=100,
                     qpCliqueOutput=NULL, density.digits=0,
                     logscale.clqsize=FALSE,
                     titleclq="maximum clique size as function of threshold",
                     verbose=FALSE) {
  n.var <- nrow(nrrMatrix)
  n.adj <- n.var*(n.var-1)/2

  if (is.null(qpCliqueOutput)) {

    if (length(breaks) == 1) {
      len <- threshold.lim[2] - threshold.lim[1]
      br <- seq(threshold.lim[1],threshold.lim[2],by=len/breaks)
    } else {
      br <- breaks
      threshold.lim <- range(br)
    }

    maxclqszeunderN <- 0
    thrmaxclqszeunderN <- 0
    mpctedclqsze <- matrix(rep(0,length(br)*3),nrow=length(br),ncol=3)
    colnames(mpctedclqsze) <- c("density","clqsize","threshold")

    ## by now we coerce this to a regular matrix
    nrrMatrix <- as(nrrMatrix, "matrix")

    for (i in 1:length(br)) {
      if (verbose) {
        cat(paste("break: ",i,sep=""))
        cat("\n")
      }
      threshold <- br[i]
      nrrMatrix[is.na(nrrMatrix)] <- 1.0 # non-available NRRs imply no edges
      A <- nrrMatrix <= threshold
      diag(A) <- FALSE # if the threshold is 1.0 the resulting qp-graph
                       # will be the complete undirected graph but
                       # still it should have no loops
      n.edg <- sum(A) / 2
      dimnames(A) <- list(1:length(A[,1]), 1:length(A[1,]))
      maxsize <- qpCliqueNumber(A, exact.calculation, approx.iter=approx.iter,
                                  verbose=verbose)
      mpctedclqsze[i,] <- c((n.edg / n.adj) * 100, maxsize, threshold)
      if (!is.na(n)) {
        if (maxsize > maxclqszeunderN && maxsize < n) {
          maxclqszeunderN <- maxsize
          thrmaxclqszeunderN <- threshold
        }
      }
    }

    if (is.na(n))
      maxclqszeunderN <- thrmaxclqszunderN <- NA

  } else {
    br <- qpCliqueOutput$data[,3]
    mpctedclqsze <- qpCliqueOutput$data
    thrmaxclqszeunderN <- qpCliqueOutput$threshold
    maxclqszeunderN <- qpCliqueOutput$clqsizeunderN
    n <- qpCliqueOutput$N
    exact.calculation <- qpCliqueOutput$exact.calculation 
  }

  linetype <- 1
  label <- "exact maximum clique size"
  if (exact.calculation == FALSE) {
    linetype <- 2
    label <- "lower bound on the maximum clique size"
  }

  if (plot == TRUE) {
    logscale <- ""
    if (logscale.clqsize == TRUE) {
      logscale <- "y"
    }
    plot(br,mpctedclqsze[,2],type="o",xlim=threshold.lim,
         ylim=range(1,n,mpctedclqsze[,2],na.rm=TRUE),lwd=2,lty=linetype,
         xlab="threshold",ylab="maximum clique size",main=titleclq,log=logscale)
    if (!is.na(n))
      abline(h=n,col="red",lwd=2)
    pct <- round(mpctedclqsze[,1],digits=density.digits)
    text(br,mpctedclqsze[,2],lab=paste(pct,"%",sep=""),pos=1,cex=.7)
    legend(min(threshold.lim),max(n,mpctedclqsze[,2],na.rm=TRUE),label,lty=linetype,lwd=2,pch=1)
  }

  m <- mpctedclqsze[,c(3,2)]
  m[,2] <- m[,2] / n.var
  f <- approxfun(m)
  area <- integrate(f,min(m[,1]),max(m[,1]))

  invisible(list(data=mpctedclqsze,complexity=area$value,threshold=thrmaxclqszeunderN,
                 clqsizeunderN=maxclqszeunderN,n=n,exact.calculation=exact.calculation))
}



## function: qpBoundary
## purpose: calculate and plot the maximum clique size as function of the
##          non-rejection rate
## parameters: nrrMatrix - matrix of non-rejection rates
##             n - sample size
##             threshold.lim - range of threshold values
##             breaks - either a number of threshold bins or a vector of
##                      threshold breakpoints
##             vertexSubset - subset of vertices for which we calculate
##                            their maximum boundary with respect to all other vertices
##             plot - flag that when set it makes a plot of the result
##             qpBoundaryOutput - output from a previous call to qpClique, this
##                              allows one to plot the result changing some of
##                              the plotting parameters without having to do
##                              the calculation again
##             density.digits - number of digits in the reported graph densities
##             logscale.bdsize - whether the boundary size in the y-axis
##                               should be plotted in logarithmic scale
##             titlebd - title to be shown in the plot
##             verbose - show progress when doing the calculation
## return: a list with the maximum boundary and the graph density as function
##         of the threshold, an estimate of the complexity of the resulting
##         qp-graphs across the thresholds, the threshold on the non-rejection
##         rate that provides a maximum clique size strictly smaller than the
##         sample size N, and the resulting maximum clique size

qpBoundary <- function(nrrMatrix, n=NA, threshold.lim=c(0,1), breaks=5, vertexSubset=NULL,
                       plot=TRUE, qpBoundaryOutput=NULL, density.digits=0,
                       logscale.bdsize=FALSE,
                       titlebd="Maximum boundary size as function of threshold",
                       verbose=FALSE) {
  n.var <- nrow(nrrMatrix)
  n.adj <- n.var*(n.var-1)/2

  if (is.null(qpBoundaryOutput)) {

    if (length(breaks) == 1) {
      len <- threshold.lim[2] - threshold.lim[1]
      br <- seq(threshold.lim[1],threshold.lim[2],by=len/breaks)
    } else {
      br <- breaks
      threshold.lim <- range(br)
    }

    maxbdszeunderN <- 0
    thrmaxbdszeunderN <- 0
    mpctedbdsze <- matrix(rep(0,length(br)*3),nrow=length(br),ncol=3)
    colnames(mpctedbdsze) <- c("density","bdsize","threshold")

    ## by now we coerce this to a regular matrix
    nrrMatrix <- as(nrrMatrix, "matrix")

    for (i in 1:length(br)) {
      if (verbose) {
        cat(paste("break: ",i,sep=""))
        cat("\n")
      }
      threshold <- br[i]
      nrrMatrix[is.na(nrrMatrix)] <- 1.0 # non-available NRRs imply no edges
      A <- nrrMatrix <= threshold
      diag(A) <- FALSE # if the threshold is 1.0 the resulting qp-graph
                       # will be the complete undirected graph but
                       # still it should have no loops
      n.edg <- sum(A) / 2
      ## dimnames(A) <- list(1:length(A[,1]), 1:length(A[1,])) WHY THIS ???
      maxsize <- max(rowSums(A))
      if (!is.null(vertexSubset))
        maxsize <- max(rowSums(A[vertexSubset, ]))

      mpctedbdsze[i,] <- c((n.edg / n.adj) * 100, maxsize, threshold)
      if (!is.na(n)) {
        if (maxsize > maxbdszeunderN && maxsize < n) {
          maxbdszeunderN <- maxsize
          thrmaxbdszeunderN <- threshold
        }
      }
    }

    if (is.na(n))
      maxbdszeunderN <- thrmaxbdszunderN <- NA

  } else {
    br <- qpBoundaryOutput$data[,3]
    mpctedbdsze <- qpBoundaryOutput$data
    thrmaxbdszeunderN <- qpBoundaryOutput$threshold
    maxbdszeunderN <- qpBoundaryOutput$bdsizeunderN
    n <- qpBoundaryOutput$n
    vertexSubset <- qpBoundaryOutput$vertexSubset
  }

  linetype <- 1
  label <- "Maximum boundary size"

  if (plot == TRUE) {
    logscale <- ""
    if (logscale.bdsize == TRUE) {
      logscale <- "y"
    }
    plot(br,mpctedbdsze[,2],type="o",xlim=threshold.lim,
         ylim=range(1,n,mpctedbdsze[,2],na.rm=TRUE),lwd=2,lty=linetype,
         xlab="threshold",ylab="Maximum boundary size",main=titlebd,log=logscale)
    if (!is.na(n))
      abline(h=n,col="red",lwd=2)
    pct <- round(mpctedbdsze[,1],digits=density.digits)
    text(br,mpctedbdsze[,2],lab=paste(pct,"%",sep=""),pos=1,cex=.7)
    legend(min(threshold.lim),max(n,mpctedbdsze[,2],na.rm=TRUE),label,lty=linetype,lwd=2,pch=1)
  }

  m <- mpctedbdsze[,c(3,2)]
  m[,2] <- m[,2] / n.var
  f <- approxfun(m)
  area <- integrate(f,min(m[,1]),max(m[,1]))

  invisible(list(data=mpctedbdsze, complexity=area$value, threshold=thrmaxbdszeunderN,
                 bdsizeunderN=maxbdszeunderN, n=n, vertexSubset=vertexSubset))
}



## function: qpGetCliques
## purpose: find the set of (maximal) cliques of an undirected graph.
##          it calls the C compiled faster code of the functions
##          from the cliquer library implementing the algorithms from
##
##          Ostergard, PRJ. A fast algorithm for the maximum clique problem
##          Discrete Appl. Math., 120:195-205, 2002
##          http://users.tkk.fi/~pat/cliquer.html
##
## parameters: g - either a graphNEL object or an adjacency matrix of the graph
##             clqspervtx - store the indices of the cliques where each vertex
##                          belongs to in the first p entries (|V|=p) of the
##                          returned list
##             verbose - show progress on the clique search
## return: a list of maximal cliques

qpGetCliques <- function(g, clqspervtx=FALSE, verbose=TRUE) {

  if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    ## require(graph)
    if (graph::edgemode(g) != "undirected")
      stop("g should be an undirected graph\n")

    A <- as(g, "matrix") == 1
  } else if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    A <- g
    p <- (d <- dim(A))[1]
    if (p != d[2])
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!isSymmetric(A))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g should be either a graphNEL object, a graphAM object, a graphBAM object or a boolean adjacency matrix\n")

  p <- dim(A)[1]
  nEd <- sum(A)/2
  if (nEd == 0) {
    clqlst <- as.list(1:p)
    if (clqspervtx)
      clqlst <- c(clqlst, as.list(1:p))
    return(clqlst)
  }

  if (nEd == (p*(p-1))/2) {
    clqlst <- list(1:p)
    if (clqspervtx)
      clqlst <- c(clqlst, as.list(rep(1, times=p)))
    return(list(1:p))
  }

  return(qpgraph:::.qpFastCliquerGetCliques(as.matrix(A), clqspervtx=clqspervtx,
                                            verbose=verbose))
}



## function: qpUpdateCliquesRemoving
## purpose: update the set of maximal cliques of an undirected graph after
##          removing a given edge in order to obtain the corresponding set of
##          maximal cliques after the edge removal. it calls the C compiled faster
##          code of a function implementing the algorithm from
##
##          Stix, V. Finding all maximal cliques in dynamic graphs
##          Comput. Optimization and Appl., 27:173-186, 2004.
##
## parameters: g - either a graphNEL object or an adjacency matrix of the graph
##             clqslst - list of cliques of the graph encoded in g. this list should
##                       start on element n+1 (for n vertices) while between elements
##                       1 to n there should be references to the cliques to which each
##                       of the 1 to n vertices belong to (i.e., the output of
##                       qpGetCliques() with parameter clqspervtx=TRUE 
##             v - vertex of the edge being removed
##             w - vertex of the edge being removed
##             verbose - show progress on the clique search
## return: an updated list of maximal cliques

qpUpdateCliquesRemoving <- function(g, clqlst, v, w, verbose=TRUE) {
  if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    ## require(graph)
    if (graph::edgemode(g) != "undirected")
      stop("g should be an undirected graph\n")

    A <- as(g, "matrix") == 1
  } else if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    A <- g
    p <- (d <- dim(A))[1]
    if (p != d[2])
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!isSymmetric(A))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g should be either a graphNEL object, a graphAM object, a graphBAM object or a boolean adjacency matrix\n")

  if (is.character(v)) {
    v <- match(v, colnames(A))
    if (is.na(v))
      stop("vertex ", v, " does not match any vertex in 'g'\n")
  }

  if (is.character(w)) {
    w <- match(w, colnames(A))
    if (is.na(w))
      stop("vertex ", w, " does not match any vertex in 'g'\n")
  }

  return(qpgraph:::.qpFastUpdateCliquesRemoving(as.matrix(A), clqlst, v, w, verbose=verbose))
}



## function: qpIPF
## purpose: implement the Iterative Proportional Fitting (IPF) algorithm to
##          perform maximum likelihood estimation of the sample covariance
##          matrix given the independence constraints from an input list of
##          (maximal) cliques.
##          Part of the R code below has been borrowed from an implementation
##          by Graham Wills in June of 1992
## parameters: vv - input matrix (usually the sample variance-covariance matrix)
##             clqlst - list of (maximal) cliques
##             tol - tolerance under which the main loop stops
##             verbose - when set to TRUE the algorithm shows the successive
##                       precision achieved at each iteration
##             R.code.only - flag set to FALSE when using the C implementation
## return: the input matrix adjusted to the constraints of the list of cliques

qpIPF <- function(vv, clqlst, tol = 0.001, verbose = FALSE,
                  R.code.only = FALSE) {

  ## this should be changed so that the rest of the algorithm deals with dspMatrix matrices
  vv <- as(vv, "matrix")

  if (!R.code.only) {
    return(qpgraph:::.qpFastIPF(vv, clqlst, tol, verbose))
  }

  if (verbose) {
    n.var <- nrow(vv)
    if (clqlst[[1]][1] <= n.var) {
      n.var <- 0
    }
    cat(paste(paste("qpIPF: ",length(clqlst)-n.var),"cliques\n"))
  }

  V <- diag(length(vv[, 1]))
  precision <- 1
  while(precision > tol) {
    Vold <- V
    V <- qpgraph:::.qpIPFpass(vv, V, clqlst)
    precision <- max(abs(V - Vold))
    if (verbose)
      cat("qpIPF: precision =", precision, "\n")
  }

  return(as(V, "dspMatrix"))
}



## function: qpHTF
## purpose: implement the algorithm of Hastie, Tibshirani and Friedman to
##          perform maximum likelihood estimation of the sample covariance
##          matrix given the independence constraints from an input
##          undirected graph
##          Part of the R code below has been borrowed from an implementation
##          by Giovanni Machetti in the 'ggm' package (thanks Giovanni!!)
## parameters: S - sample variance-covariance matrix
##             g - input undirected graph
##             tol - tolerance under which the main loop stops
##             verbose - when set to TRUE the algorithm shows the successive
##                       precision achieved at each iteration
##             R.code.only - flag set to FALSE when using the C implementation
## return: the input matrix adjusted to the constraints of the list of cliques

qpHTF <- function(S, g, tol = 0.001, verbose = FALSE,
                  R.code.only = FALSE) {

  ## this should be changed so that the rest of the algorithm deals with dspMatrix matrices
  S <- as(S, "matrix")

  A <- NULL
  n.var <- NULL
  var.names <- NULL
  if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    n.var <- nrow(g)
    var.names <- rownames(g)
    if (is.null(rownames(var.names)))
      var.names <- 1:n.var
    A <- g
    if (class(A[1, 1]) != "logical")
      A <- A == 1 ## get a logical adjacency matrix
    ## by now we have to coerce the adjacency matrix to a regular matrix
    ## but in the future this should be working with the more memory-efficient dspMatrix class
    A <- as.matrix(A)
  }
  if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    n.var <- length(graph::nodes(g))
    var.names <- nodes(g)
    A <- as(g, "matrix") == 1 ## get a logical adjacency matrix
  }

  if (is.null(n.var))
    stop("'g' is neither an adjacency matrix, a graphNEL, a graphAM or a graphBAM object.\n")

  if (verbose)
    cat("qpHTF: maximum boundary =", max(rowSums(A)), "\n") 


  if (!R.code.only) {
    return(qpgraph:::.qpFastHTF(S, A, tol, verbose))
  }

  ppct <- -1
  pb <- NULL
  if (verbose)
    pb <- txtProgressBar(style=3)

  W <- Wprev <- S
  precision <- 1
  while (precision > tol) {
    Wprev <- W
    for (i in 1:n.var) {
      W11 <- W[-i, -i, drop=FALSE]
      s12 <- S[-i, i, drop=FALSE]
      Ai <- A[i, ]
      Ai <- Ai[-i]
      beta <- rep(0, n.var-1)
      beta[Ai] <- solve(W11[Ai, Ai, drop=FALSE], s12[Ai, , drop=FALSE])
      w <- W11 %*% beta
      W[-i, i] <- W[i, -i] <- w

      pct <- floor((i * 100) / n.var)
      if (pct != ppct && verbose) {
        setTxtProgressBar(pb, pct/100)
        ppct <- pct
      }
    }
    precision <- max(abs(W - Wprev))
    if (verbose)
      cat("\nqpHTF: precision =", precision, "\n")
  }

  return(as(W, "dspMatrix"))
}



## function: qpPAC
## purpose: for a given undirected graph in an adjacency matrix estimate the
##          partial correlation coefficient (PAC) and its corresponding p-value
##          for each edge in the graph
## parameters: X - data set from where to estimate the PACs
##             g - either a graphNEL object or an adjacency matrix of the graph
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             return.K - flag that when set to TRUE the function also returns
##                        the concentration matrix K; if FALSE (default) does not
##                        return K
##             tol - maximum tolerance in the application of the IPF algorithm
##             matrix.completion - algorithm to perform matrix completion operations
##             verbose - flag that when set to TRUE the IPF algorithm
##                       shows the convergence progression
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two matrices, one with the estimates of the PACs and
##         the other with their p-values

setGeneric("qpPAC", function(X, ...) standardGeneric("qpPAC"))

# X comes as an ExpressionSet object
setMethod("qpPAC", signature(X="ExpressionSet"),
          function(X, g, return.K=FALSE, tol=0.001, matrix.completion=c("HTF", "IPF"),
                   verbose=TRUE, R.code.only=FALSE) {
            X <- t(Biobase::exprs(X))
            qpgraph:::.qpPAC(X, g, return.K, tol, matrix.completion, verbose, R.code.only)
          })

# X comes as a data frame
setMethod("qpPAC", signature(X="data.frame"),
          function(X, g, return.K=FALSE, long.dim.are.variables=TRUE,
                   tol=0.001, matrix.completion=c("HTF", "IPF"), verbose=TRUE,
                   R.code.only=FALSE) {
            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)
            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpPAC(X, g, return.K, tol, matrix.completion, verbose, R.code.only)
          })

          
# data comes as a matrix
setMethod("qpPAC", signature(X="matrix"),
          function(X, g, return.K=FALSE, long.dim.are.variables=TRUE,
                   tol=0.001, matrix.completion=c("HTF", "IPF"), verbose=TRUE,
                   R.code.only=FALSE) {
            if (long.dim.are.variables &&
              sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpPAC(X, g, return.K, tol, matrix.completion, verbose, R.code.only)
          })

.qpPAC <- function(X, g, return.K=FALSE, tol=0.001, matrix.completion=c("HTF", "IPF"),
                   verbose=TRUE, R.code.only=FALSE) {
  matrix.completion <- match.arg(matrix.completion)

  A <- matrix(FALSE, nrow=ncol(X), ncol=ncol(X), dimnames=list(colnames(X), colnames(X)))
  if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    ## require(graph)
    if (graph::edgemode(g) != "undirected")
      stop("g should be an undirected graph\n")

    Ag <- as(g, "matrix") == 1
    if (any(is.na(match(rownames(Ag), colnames(X)))))
      stop("some variables in the graph 'g' do not match the variables in the data")

    A[rownames(Ag), colnames(Ag)] <- Ag
  } else if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    A <- g
    p <- (d <- dim(A))[1]
    if (p != d[2])
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!isSymmetric(A))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g should be either a graphNEL object or a boolean adjacency matrix\n")

  var.names <- colnames(X)
  n.var <- ncol(X)
  N <- nrow(X)

  ## calculate sample covariance matrix
  S <- qpCov(X)

  ## ensure rows and columns follow the same order
  A <- A[rownames(S), colnames(S)]

  if (matrix.completion == "IPF") {
    ## get the list of maximal cliques
    clqlst <- qpGetCliques(A, verbose=verbose)

    ## get a maximum likelihood estimate of the sample covariance matrix
    ## using the clique list and the IPF algorithm
    Shat <- qpIPF(S, clqlst, tol=tol, verbose=verbose, R.code.only=R.code.only)
  } else
    Shat <- qpHTF(S, A, tol=tol, verbose=verbose, R.code.only=R.code.only)

  ## estimate partial correlation coefficients and their standard errors

  K <- solve(Shat)
  SE <- qpgraph:::.qpEdgePACSE(Shat, A, R.code.only=R.code.only)

  ## return matrices of partial correlations, standard errors
  ## and p-values for every edge

  C <- N * (K^2 / SE)
  rho_coef <- qpK2ParCor(K)
  p.values <- 1 - pchisq(C, df=1)
  dimnames(rho_coef) <- dimnames(p.values) <- list(var.names, var.names)

  list2return <- list()

  if (return.K)
    list2return <- list(R=as(forceSymmetric(rho_coef), "dspMatrix"),
                        P=as(forceSymmetric(p.values), "dspMatrix"), K=Matrix(K))
  else
    list2return <- list(R=as(forceSymmetric(rho_coef), "dspMatrix"),
                        P=as(forceSymmetric(p.values), "dspMatrix"))

  return(list2return)
}



## function: qpPCC
## purpose: estimate pairwise Pearson correlation coefficients (PCCs) between all
##         pairs of variables
## parameters: X - data set from where to estimate the PCCs
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
## return: a list with two matrices, one with the estimated PCCs and the
##         other with their p-values

setGeneric("qpPCC", function(X, ...) standardGeneric("qpPCC"))

# X comes as an ExpressionSet object
setMethod("qpPCC", signature(X="ExpressionSet"),
          function(X) {
            X <- t(Biobase::exprs(X))
            qpgraph:::.qpPCC(X)
          })

# X comes as a data frame
setMethod("qpPCC", signature(X="data.frame"),
          function(X, long.dim.are.variables=TRUE) {
            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)
            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpPCC(X)
          })

          
# X comes as a matrix
setMethod("qpPCC", signature(X="matrix"),
          function(X, long.dim.are.variables=TRUE) {
            if (long.dim.are.variables &&
              sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)
            qpgraph:::.qpPCC(X)
          })

.qpPCC <- function(X) {

  var.names <- colnames(X)
  N <- nrow(X)

  ## calculate sample covariance matrix
  S <- qpCov(X)

  ## estimate PCCs by scaling the covariance matrix
  ## somehow Matrix::cov2cor() refuses to scale non-positive definite matrices stored as dspMatrix objects (?)
  R <- as(Matrix::cov2cor(as.matrix(S)), "dspMatrix")

  ## calculate t-statistics
  T <- (N - 2) / (1 - R*R)
  diag(T) <- (N - 2) * 100000 # just to get 0 p-values on the diagonal
  T <- R * sqrt(T)

  ## calculate two-sided p-values
  p <- pt(as.matrix(T), df=N - 2)
  P <- as(2 * pmin(p, 1 - p), "dspMatrix")

  list(R=R, P=P)
}






## function: qpK2ParCor
## purpose: obtain the partial correlation coefficients from a given
##          concentration matrix
## parameters: K - concentration matrix
## return: a matrix with the partial correlation coefficients

qpK2ParCor <- function(K) {
  R <- -cov2cor(K)
  diag(R) <- 1
  return(R)
}



## function: qpPrecisionRecall
## purpose: calculate the precision-recall curve for a given measure of
##          association with respect to some given reference graph
## parameters: measurementsMatrix - matrix containing the measure of association
##                                  between all pairs of variables
##             refGraph - a reference graph from which to calculate the precision-recall
##                        curve provided either as an adjacency matrix, a two-column matrix
##                        of edges, a graphNEL object, a graphAM object or a graphBAM object
##             decreasing - logical; if TRUE then the measurements are ordered
##                          in decreasing order; if FALSE then in increasing
##                          order
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             recallSteps - steps of the recall on which to calculate precision
## return: a matrix where rows correspond to recall steps and columns correspond,
##         respetively, to the actual recall, the precision, the number of true
##         positives at that recall rate and the threshold score that yields that
##         recall rate

qpPrecisionRecall <- function(measurementsMatrix, refGraph, decreasing=TRUE,
                              pairup.i=NULL, pairup.j=NULL,
                              recallSteps=seq(0.0, 1.0, by=0.1)) {

  if (class(measurementsMatrix) != "matrix" && class(measurementsMatrix) != "dspMatrix" &&
      class(measurementsMatrix) != "dgeMatrix")
    stop("measurementsMatrix should be a numerical matrix\n")

  p <- (d <- dim(measurementsMatrix))[1]
  if (p != d[2])
    stop("measurementsMatrix should be a squared matrix\n")

  if (class(refGraph) != "data.frame" && class(refGraph) != "matrix" &&
      class(refGraph) != "graphNEL" && class(refGraph) != "graphAM" &&
      class(refGraph) != "graphBAM" &&
      length(grep("Matrix", class(refGraph))) == 0)
    stop("refGraph should be provided either as an adjacency matrix, a two-column matrix of edges, a graphNEL object, a graphAM object or a graphBAM object\n")

  if (class(refGraph) == "data.frame" || class(refGraph) == "matrix" ||
      length(grep("Matrix", class(refGraph))) > 0) {
    p <- (d <- dim(refGraph))[1]
    if (p != d[2] && ncol(refGraph) != 2)
      stop("If refGraph is a matrix then it should be either a squared adjacency matrix or a two-column matrix with rows corresponding to edges \n")

    if (p != d[2] && ncol(refGraph) == 2) {
      if (class(refGraph[1, 1]) == "character") {
        refGraph <- cbind(match(refGraph[, 1], rownames(measurementsMatrix)),
                          match(refGraph[, 2], rownames(measurementsMatrix)))
        if (any(is.na(refGraph)))
          stop("Some of the identifiers in refGraph do not correspond to row and column names in measurementsMatrix\n")
      }

      refA <- matrix(FALSE, nrow=nrow(measurementsMatrix),
                            ncol=ncol(measurementsMatrix))
      refA[as.matrix(refGraph)] <- TRUE
      refA <- refA | t(refA)
      rownames(refA) <- colnames(refA) <- rownames(measurementsMatrix)
    } else {
      refA <- as.matrix(refGraph)
      if (nrow(measurementsMatrix) != nrow(refA) ||
          ncol(measurementsMatrix) != ncol(refA))
        stop("measurementsMatrix and refGraph should have the same dimensions\n")
    }

  } else { ## graphNEL, graphAM or graphBAM
    if (any(is.na(match(graph::nodes(refGraph), rownames(measurementsMatrix)))) ||
        length(graph::nodes(refGraph) != dim(measurementsMatrix)[1]))
      stop("The vertex set in refGraph does not correspond to the row and column names in measurementsMatrix\n")
    refA <- as(refGraph, "matrix") == 1
    if (!all(match(rownames(refA), rownames(measurementsMatrix)) == 1:dim(refA)[1])) { ## force matching the vertex order
      refA <- refA[match(rownames(measurementsMatrix), rownames(refA)), ] ## re-order rows
      refA <- refA[, match(colnames(measurementsMatrix), colnames(refA)), ] ## re-order columns
    }
  }

  n.var <- nrow(measurementsMatrix)

  ## by now we have to coerce it to a regular matrix
  ## but hopefully in the near future we can do [<- as well on a dspMatrix
  measurementsMatrix <- as.matrix(measurementsMatrix)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(measurementsMatrix)))
      stop("when using pairup.i and pairup.j, measurementsMatrix should have row and column names\n")

    var.names <- colnames(measurementsMatrix)
    pairup.i <- match(pairup.i, var.names)
    if (sum(is.na(pairup.i)) > 0)
      stop("pairup.i is not a subset of the variables forming the data\n")
    pairup.j <- match(pairup.j, var.names)
    if (sum(is.na(pairup.j)) > 0)
      stop("pairup.j is not a subset of the variables forming the data\n")

    pairup.ij.int <- intersect(pairup.i, pairup.j)
    pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
    pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

    nomeasurementsMask <- matrix(FALSE, nrow=n.var, ncol=n.var)
    nomeasurementsMask[as.matrix(
                       expand.grid(pairup.ij.int,
                                   c(pairup.i.noint, pairup.j.noint)))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.i.noint, pairup.j.noint))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.ij.int, pairup.ij.int))] <- TRUE
    diag(nomeasurementsMask) <- FALSE
    nomeasurementsMask <- nomeasurementsMask | t(nomeasurementsMask)
    nomeasurementsMask <- !nomeasurementsMask
    measurementsMatrix[nomeasurementsMask] <- NA
  }

  upperTriRow <- row(measurementsMatrix)[upper.tri(measurementsMatrix) &
                                         !is.na(measurementsMatrix)]
  upperTriCol <- col(measurementsMatrix)[upper.tri(measurementsMatrix) &
                                         !is.na(measurementsMatrix)]
  measurementsUpperTriMatrix <- measurementsMatrix[upper.tri(measurementsMatrix) &
                                                   !is.na(measurementsMatrix)]
  orderedMeasurements <- sort(measurementsUpperTriMatrix,index.return=TRUE,
                              decreasing=decreasing)
  edgeRnk <- cbind(upperTriRow[orderedMeasurements$ix],
                   upperTriCol[orderedMeasurements$ix],
                   orderedMeasurements$x)

  lenRnk <- dim(edgeRnk)[1]
  total_positives <- sum(refA[upper.tri(refA) & !is.na(measurementsMatrix)])

  status <- refA[as.matrix(edgeRnk[,c(1,2)])]
  status_tp <- rep(0, length(status))
  status_tp[status] <- 1:total_positives
  status_tp[which.max(status_tp):length(status_tp)] <- total_positives
  preRec <- matrix(0, nrow=length(recallSteps), ncol=5)
  colnames(preRec) <- c("Recall","Precision", "TP", "FP", "ScoreThreshold")

  for (i in 1:length(recallSteps)) {
    tp <- round(recallSteps[i] * total_positives, digits=0)
    rnkPos <- max(c(0, (1:lenRnk)[status_tp <= tp & status_tp != 0]))
    fp <- rnkPos - tp

    actualRecall <- tp / total_positives
    if (tp + fp > 0)
      precision <- tp / (tp + fp)
    else
      precision <- 1.0

    if (rnkPos > 0)
      scoreThreshold <- edgeRnk[rnkPos, 3]
    else
      scoreThreshold <- NA

    preRec[i, ] <- c(actualRecall, precision, tp, fp, scoreThreshold)
  }

  return(preRec)
}



## function: qpPRscoreThreshold
## purpose: calculate the score threshold at a given precision or recall
##          level from a given precision-recall curve
## parameters: preRecFun - precision-recall function (output from
##                         qpPrecisionRecall)
##             level - recall or precision level
##             recall.level - logical; if TRUE then it is assumed that the value
##                            given in the level parameter corresponds to
##                            a desired level of recall; if FALSE then it
##                            is assumed a desired level of precision
##             max.score - maximum score given by the method that produced
##                         the precision-recall function to an association
## return: the score threshold at which a given level of precision or recall
##         is attained by the given precision-recall function

qpPRscoreThreshold <- function(preRecFun, level, recall.level=TRUE,
                               max.score=9999999) {
  preRecFun[1, "ScoreThreshold"] <- max.score

  levelCol <- "Recall"
  i1 <- max((1:length(preRecFun[,1]))[preRecFun[,"Recall"] <= level])

  if (!recall.level) {
    i1 <- max((1:length(preRecFun[,1]))[preRecFun[,"Precision"] >= level])
    levelCol <- "Precision"
  }

  i0 <- i1 + 1
  x1 <- preRecFun[i1, levelCol]
  x0 <- preRecFun[i0, levelCol]
  y1 <- preRecFun[i1, "ScoreThreshold"]
  y0 <- preRecFun[i0, "ScoreThreshold"]

  y <- y0 + (level - x0) * (y1 - y0) / (x1 - x0)

  return(y)
}



## function: qpImportNrr
## purpose: import non-rejection rates from a flat file
##
## parameters: filename - name of the flat file with the non-rejection rates
##             nTests - number of tests performed in the estimation of those
##                      non-rejection rates
## return: a matrix of non-rejection rates

qpImportNrr <- function(filename, nTests) {
  nrr <- as.matrix(read.table(filename))

  n.var <- max(nrr[,c(1,2)]) + 1
  nrrMatrix <- matrix(as.double(0), nrow=n.var, ncol=n.var)
  nrrMatrix[nrr[,c(1,2)]+1] <- nrr[,3] / nTests
  nrrMatrix <- nrrMatrix + t(nrrMatrix)
  diag(nrrMatrix) <- NA

  return(nrrMatrix)
}



## function: qpFunctionalCoherence
## purpose: estimate functional coherence of a transcripcional regulatory network
##          represented by means of an undirected graph encoded by an adjacency
##          matrix and of a set of transcription factor genes. In these
##          calculations Gene Ontology (GO) annotations are employed through a
##          given annotation .db package for the Entrez Gene IDs associated to
##          the rows and columns of the adjacency matrix.
## parameters: object - adjacency matrix of the undirected graph representing the
##                 transcriptional regulatory network
##             TFgenes - vector of transcription factor gene names (matching the
##                       genes at the rows and column names of A)
##             geneUniverse - vector of all genes relevant to the analysis
##             chip - name of the .db package containing the GO annotations
##             minRMsize - minimum size of the target gene set in each regulatory
##                         module where functional enrichment will be calculated
##                         and thus where functional coherence will be estimated
##             verbose - logical; if TRUE the function will show progress on the
##                       calculations; if FALSE will remain quiet (default)
##             clusterSize - size of the cluster of processors to do calculations
##                           in parallel via 'snow' and 'rlecuyer'
## return: a list with three slots, a first one containing the transcriptional
##         regulatory network as a list of regulatory modules and their targets,
##         a second one containing this same network but including only those
##         modules with GO BP annotations and a third one consisting of a vector
##         of functional coherence values

setGeneric("qpFunctionalCoherence", function(object, ...) standardGeneric("qpFunctionalCoherence"))

## the input object is a lsCMatrix adjacency matrix
setMethod("qpFunctionalCoherence",
          signature(object="lsCMatrix"),
          function(object, TFgenes, geneUniverse=rownames(object), chip, minRMsize=5, removeGOterm="transcription", verbose=FALSE, clusterSize=1)
            qpFunctionalCoherence(as(object, "matrix"), TFgenes, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize))

## the input object is a lspMatrix adjacency matrix
setMethod("qpFunctionalCoherence",
          signature(object="lspMatrix"),
          function(object, TFgenes, geneUniverse=rownames(object), chip, minRMsize=5, removeGOterm="transcription", verbose=FALSE, clusterSize=1)
            qpFunctionalCoherence(as(object, "matrix"), TFgenes, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize))

## the input object is a lsyMatrix adjacency matrix
setMethod("qpFunctionalCoherence",
          signature(object="lsyMatrix"),
          function(object, TFgenes, geneUniverse=rownames(object), chip, minRMsize=5, removeGOterm="transcription", verbose=FALSE, clusterSize=1)
            qpFunctionalCoherence(as(object, "matrix"), TFgenes, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize))

## the input object is a regular adjacency matrix
setMethod("qpFunctionalCoherence",
          signature(object="matrix"),
          function(object, TFgenes, geneUniverse=rownames(object), chip, minRMsize=5, removeGOterm="transcription", verbose=FALSE, clusterSize=1) {

  if (!qpgraph:::.qpIsPackageLoaded("GOstats"))
    stop("qpFunctionalCoherence() requires first loading Bioconductor package GOstats")

  if (clusterSize > 1 &&
      (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
    stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

  if (is.null(colnames(object)) || is.null(rownames(object)))
    stop("the adjacency matrix contained in the 'object' argument should have row and column names corresponding to the gene IDs")

  if (class(object[1,1]) != "logical" && class(object[1,1]) != "numeric" && class(object[1,1]) != "integer")
    stop("the adjacency matrix should be either logical or binary")

  if (class(object[1,1]) == "numeric" || class(object[1,1]) == "integer")
    object <- object == 1

  if (length(TFgenes) < 1)
    stop("TFgenes should contain at least one transcription factor gene\n")

  if (!is.character(TFgenes))
    stop("gene identifiers in TFgenes should belong to the class character\n")

  if (sum(is.na(match(TFgenes, geneUniverse))) > 0)
    stop("TFgenes is not a subset from the genes comprising the gene universe\n")

  TFgenes_i <- match(TFgenes, geneUniverse)
  txRegNet <- lapply(TFgenes_i, function(x) geneUniverse[object[as.integer(x), ]])
  names(txRegNet) <- TFgenes

  return(qpgraph:::.qpFunctionalCoherence(txRegNet, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize))
})

## the input object is a list of regulatory modules
setMethod("qpFunctionalCoherence",
          signature(object="list"),
          function(object, geneUniverse=unique(c(names(object), unlist(object, use.names=FALSE))),
                   chip, minRMsize=5, removeGOterm="transcription", verbose=FALSE, clusterSize=1) {

  if (clusterSize > 1 &&
      (!qpgraph:::.qpIsPackageLoaded("rlecuyer") || !qpgraph:::.qpIsPackageLoaded("snow")))
    stop("Using a cluster (clusterSize > 1) requires first loading packages 'snow' and 'rlecuyer'\n")

  TFgenes <- names(object)

  if (length(TFgenes) < 1)
    stop("names of the entries in the input list should correspond to the transcription factor gene identifiers\n")

  if (sum(is.na(match(TFgenes, geneUniverse))) > 0)
    stop("TFgenes is not a subset from the genes comprising the gene universe\n")

  return(qpgraph:::.qpFunctionalCoherence(object, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize))
})

.qpFunctionalCoherence <- function(txRegNet, geneUniverse, chip, minRMsize, removeGOterm, verbose, clusterSize) {
  TFgenes <- names(txRegNet)
  regModuleSize <- unlist(lapply(txRegNet, length))

  hpTest <- get("hyperGTest", mode="function")
  goGraph <- get("GOGraph", mode="function")
  simui <- get("simUI", mode="function")
  sigcat <- get("sigCategories", mode="function")

  geneBPuniverse <- qpgraph:::.qpFilterByGO(geneUniverse, chip, "BP")

  if (sum(regModuleSize >= minRMsize) == 0)
    stop(sprintf("qpFunctionalCoherence: No regulatory module has a minimum size of %d\n", minRMsize))

  if (verbose && minRMsize > 1)
    cat(sprintf("qpFunctionalCoherence: calculating GO enrichment in %d regulatory modules\n",
        length(TFgenes[regModuleSize >= minRMsize])))

  ## WARNING: THIS IS REALLY NECESSARY ONLY WHEN THE TF HAS GO ANNOTATIONS !!!!
  if (clusterSize > 1) {
    ## copying ShortRead's strategy, 'get()' are to quieten R CMD check, and for no other reason
    makeCl <- get("makeCluster", mode="function")
    clSetupRNG <- get("clusterSetupRNG", mode="function")
    clEvalQ <- get("clusterEvalQ", mode="function")
    clExport <- get("clusterExport", mode="function")
    clApply <- get("clusterApply", mode="function")
    stopCl <- get("stopCluster", mode="function")
    clCall <- get("clusterCall", mode="function")
    clOpt <- get("getClusterOption", mode="function")
    clParLapply <- get("parLapply", mode="function")

    message("Estimating functional coherence using a ", clOpt("type"),
            " cluster of ", clusterSize, " nodes\n")

    cl <- makeCl(clusterSize, type="MPI", snowlib=system.file(package="qpgraph"))
    clSetupRNG(cl)
    res <- clEvalQ(cl, require(qpgraph, quietly=TRUE))
    if (!all(unlist(res))) {
      stopCl(cl)
      stop("The package 'qpgraph' could not be loaded in some of the nodes of the cluster")
    }
    res <- clEvalQ(cl, require(GOstats, quietly=TRUE))
    if (!all(unlist(res))) {
      stopCl(cl)
      stop("The package 'GOstats' could not be loaded in some of the nodes of the cluster")
    }
    txRegNetGO <- clParLapply(cl, txRegNet[TFgenes[regModuleSize >= minRMsize]],
                              function(TFgeneTGs, geneBPuniverse, chip) {
                                TFgeneTGsWithGO <- qpgraph:::.qpFilterByGO(TFgeneTGs, chip, "BP")
                                res <- list(TGgenesWithGO=TFgeneTGsWithGO,
                                            goBPcondResult=NULL,
                                            goBPcondResultSigCat=NULL)

                                if (length(TFgeneTGsWithGO) >= minRMsize && minRMsize > 1) {
                                  goHypGparams <- new("GOHyperGParams",
                                                      geneIds=TFgeneTGsWithGO,
                                                      universeGeneIds=geneBPuniverse,
                                                      annotation=chip, ontology="BP",
                                                      pvalueCutoff=0.05, conditional=TRUE,
                                                      testDirection="over")
                                  goHypGcond <- hpTest(goHypGparams)
                                  res <- list(TGgenesWithGO=TFgeneTGsWithGO,
                                              goBPcondResult=goHypGcond,
                                              goBPcondResultSigCat=sigcat(goHypGcond))
                                }
                                res
                             }, geneBPuniverse, chip)
    txRegNetGO <- txRegNetGO[which(sapply(txRegNetGO, length) > 0)]
    stopCl(cl)
  } else {
    txRegNetGO <- list()
    for (TFgene in TFgenes[regModuleSize >= minRMsize]) {
      if (verbose)
        cat(".")
      TFgeneTGs <- txRegNet[[TFgene]]
      TFgeneTGsWithGO <- qpgraph:::.qpFilterByGO(TFgeneTGs, chip, "BP")

      txRegNetGO[[TFgene]] <- list(TGgenesWithGO=TFgeneTGsWithGO,
                                   goBPcondResult=NULL,
                                   goBPcondResultSigCat=NULL)

      if (length(TFgeneTGsWithGO) >= minRMsize && minRMsize > 1) {
        goHypGparams <- new("GOHyperGParams",
                            geneIds=TFgeneTGsWithGO,
                            universeGeneIds=geneBPuniverse,
                            annotation=chip, ontology="BP",
                            pvalueCutoff=0.05, conditional=TRUE,
                            testDirection="over")
        goHypGcond <- hpTest(goHypGparams)
        txRegNetGO[[TFgene]]$goBPcondResult <- goHypGcond
        txRegNetGO[[TFgene]]$goBPcondResultSigCat <- sigcat(goHypGcond)
      }
    }
  }

  if (verbose)
    cat(sprintf("\nqpFunctionalCoherence: calculating functional coherence in %d regulatory modules\n",
        length(names(txRegNetGO))))

  TFgenesWithGO <- qpgraph:::.qpFilterByGO(TFgenes, chip, "BP")
  TFgenesWithGO <- AnnotationDbi::mget(TFgenesWithGO, get(gsub(".db","GO",chip)))
  TFgenesWithGO <- lapply(TFgenesWithGO,
                          function(x) if (is.list(x)) {
                                        z <- sapply(x, function(x) x$Ontology);
                                        z[unique(names(z))]
                                      })
  TFgenesWithGOBP <- lapply(TFgenesWithGO,
                            function(x) if (sum(x=="BP",na.rm=TRUE) > 0) {
                                          names(x[x=="BP" & !is.na(x)])
                                        } else { NULL })

  goTermsEnv <- GOenv("TERM")
  goBPparentsEnv <- GOenv("BPPARENTS")

  TFgoTerms <- NULL
  if (!is.null(removeGOterm)) {
    ## remove from the GO annotation of the transcription factor and the target gene,
    ## if it's a single one, those GO terms matching the word in the argument
    ## 'removeGOTerm'. By default this is the word 'transcription', i.e., try to
    ## remove GO terms associated to transcriptional regulation from the transcription
    ## factor GO annotation

    goTerms <- unlist(AnnotationDbi::eapply(goTermsEnv, function(x) x@Term))
    goTermOntologies <- unlist(AnnotationDbi::eapply(goTermsEnv, function(x) x@Ontology))
    goTermBPOntology <- names(goTermOntologies[goTermOntologies == "BP"])

    TFgoTerms <- names(goTerms[grep(removeGOterm, goTerms)])
    TFgoTerms <- TFgoTerms[!is.na(match(TFgoTerms, goTermBPOntology))]
  }

  GOgraphSim <- rep(NA, length(names(txRegNetGO)))
  names(GOgraphSim) <- names(txRegNetGO)
  for (TFgene in names(txRegNetGO)) {
    if (verbose)
      cat(".")
    sUI <- NA
    # if the transcription factor has no GO Biological Process (BP) annotations
    # then the functional coherence value is NA
    if (length(TFgenesWithGOBP[[TFgene]]) > 0) {
      TFgoAnnot <- TFgenesWithGOBP[[TFgene]]
      mt <- match(TFgoAnnot, TFgoTerms)
      if (sum(!is.na(mt)) > 0)
        TFgoAnnot <- TFgoAnnot[is.na(mt)]

      # if the transcription factor has no GO BP annotations beyond
      # transcriptional regulation then the functional coherence value is NA
      if (length(TFgoAnnot) > 0) {
        txRegNetGO[[TFgene]][["TFgeneGOannot"]] <- TFgoAnnot

        # if the transcription factor has GO BP annotations beyond transcription
        # but the target gene set has no over-represented GO terms then the
        # functional coherence is 0
        if (length(txRegNetGO[[TFgene]]$goBPcondResultSigCat) == 0 && length(txRegNetGO[[TFgene]]$TGgenesWithGO) > 1)
          sUI <- 0
        else { # otherwise the functional coherence is estimated as the similarity
               # between the GO graphs associated to the transcription factor GO
               # annotations and the GO over-represented terms in the target gene set
          if (length(txRegNetGO[[TFgene]]$TGgenesWithGO) > 0) {
            gTF <- goGraph(TFgoAnnot, goBPparentsEnv)
            gTF <- removeNode("all", gTF)
            gTG <- goGraph("all", goBPparentsEnv)
            if (length(txRegNetGO[[TFgene]]$TGgenesWithGO) > 1)
              gTG <- goGraph(txRegNetGO[[TFgene]]$goBPcondResultSigCat, goBPparentsEnv)
            else { ## there's just one gene as target
              TGgene <- txRegNetGO[[TFgene]]$TGgenesWithGO
              if (length(qpgraph:::.qpFilterByGO(TGgene, chip, "BP")) > 0) {
                TGgeneWithGO <- AnnotationDbi::mget(TGgene, get(gsub(".db","GO",chip)))
                TGgeneWithGO <- lapply(TGgeneWithGO,
                                        function(x) if (is.list(x)) {
                                                      z <- sapply(x, function(x) x$Ontology);
                                                      z[unique(names(z))]
                                                    })
                TGgeneWithGOBP <- lapply(TGgeneWithGO,
                                function(x) if (sum(x=="BP",na.rm=TRUE) > 0) {
                                              names(x[x=="BP" & !is.na(x)])
                                            } else { NULL })
                TGgoAnnot <- TGgeneWithGOBP[[TGgene]]
                if (length(TGgoAnnot) > 0) {
                  mt <- match(TGgoAnnot, TFgoTerms)
                  if (sum(!is.na(mt)) > 0)
                    TGgoAnnot <- TGgoAnnot[is.na(mt)]
                  if (length(TGgoAnnot) > 0)
                    gTG <- goGraph(TGgoAnnot, goBPparentsEnv)
                  else
                    gTG <- goGraph("all", goBPparentsEnv)
                } else
                  gTG <- goGraph("all", goBPparentsEnv)
                txRegNetGO[[TFgene]][["TGgeneGOannot"]] <- TGgoAnnot
              }
            }
            gTG <- removeNode("all", gTG)

            sUI <- simui(gTF, gTG)
          }
        }
      }
    }
    GOgraphSim[TFgene] <- sUI
  }
  if (verbose)
    cat("\n")

  return(list(txRegNet=txRegNet,
              txRegNetGO=txRegNetGO,
              functionalCoherenceValues=GOgraphSim))
}



## function: qpTopPairs
## purpose: report a top number of pairs of variables from a network encoded by a given
##          graph and possibly adding annotation and other information
## parameters: measurementsMatrix - matrix of pairwise associations
##             refGraph - undirected graph of selected interactions provided either as
##                        an adjacency matrix, a graphNEL object, a graphAM object or a graphBAM object
##             n - number of pairs to report
##             file - file name to dump the pairs information as tab-separated column text
##             decreasing - logical; if TRUE then the measurements are ordered
##                          in decreasing order; if FALSE then in increasing
##                          order. It applies only when measurementsMatrix is not null
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             annotation - name of an annotation package to transform gene identifiers
##                          into gene symbols
##             fcOutput - output of qpFunctionalCoherence
##             digits - number of decimal digits reported in the association measure and
##                      functional coherence values


qpTopPairs <- function(measurementsMatrix=NULL, refGraph=NULL, n=6L, file=NULL,
                       decreasing=FALSE, pairup.i=NULL, pairup.j=NULL,
                       annotation=NULL, fcOutput=NULL, fcOutput.na.rm=FALSE,
                       digits=NULL) {

  haveMeasurements <- FALSE

  if (is.null(measurementsMatrix) && is.null(refGraph))
    stop("A proper value for either 'measurementsMatrix' or 'refGraph' should be provided\n")

  if (is.null(measurementsMatrix)) {
    if (class(refGraph) != "matrix" && class(refGraph) != "graphNEL" &&
        class(refGraph) != "graphAM" && class(refGraph) != "graphBAM" &&
        length(grep("Matrix", class(refGraph))) == 0)
      stop("refGraph should be provided either as an adjacency matrix, a graphNEL object, a graphAM object or a graphBAM object\n")

    refGraph <- as(refGraph, "matrix")

    p <- (d <- dim(refGraph))[1]
    if (p != d[2])
      stop("'measurementsMatrix' should be a squared matrix\n")

    if (is.null(rownames(refGraph)) || is.null(colnames(refGraph)))
      stop("'refGraph' should have row and column names\n")
    else {
      if (!identical(rownames(refGraph), colnames(refGraph)))
        stop("Row and column names of 'refGraph' should be the same\n")
    }

    measurementsMatrix <- matrix(0, nrow=nrow(refGraph), ncol=ncol(refGraph), dimnames=dimnames(refGraph))
  } else {
    if (class(measurementsMatrix) != "matrix" && class(measurementsMatrix) != "dspMatrix" &&
        class(measurementsMatrix) != "dgeMatrix")
      stop("'measurementsMatrix' should be a numerical matrix\n")

    p <- (d <- dim(measurementsMatrix))[1]
    if (p != d[2])
      stop("'measurementsMatrix' should be a squared matrix\n")

    if (is.null(rownames(measurementsMatrix)) || is.null(colnames(measurementsMatrix)))
      stop("'measurementsmatrix' should have row and column names\n")
    else {
      if (!identical(rownames(measurementsMatrix), colnames(measurementsMatrix)))
        stop("Row and column names of 'measurementsMatrix' should be the same\n")
    }

    haveMeasurements <- TRUE
  }

  if (is.null(refGraph)) {
    refGraph <- matrix(TRUE, nrow=nrow(measurementsMatrix),
                       ncol=ncol(measurementsMatrix),
                       dimnames=dimnames(measurementsMatrix))
  } else {
    if (class(refGraph) != "matrix" && class(refGraph) != "graphNEL" &&
        class(refGraph) != "graphAM" && class(refGraph) != "graphBAM" &&
        length(grep("Matrix", class(refGraph))) == 0)
      stop("'refGraph' should be provided either as an adjacency matrix, a graphNEL object, a graphAM object or a graphBAM object\n")

    if (class(refGraph) == "graphNEL" || class(refGraph) == "graphAM" || class(refGraph) == "graphBAM")
      refGraph <- graph::ugraph(refGraph)

    refGraph <- as(refGraph, "matrix")

    p <- (d <- dim(refGraph))[1]
    if (p != d[2])
      stop("'refGraph' should be a squared matrix\n")

    if (is.null(rownames(refGraph)) || is.null(colnames(refGraph)))
      stop("'refGraph' should have row and column names\n")
    else {
      if (!identical(rownames(refGraph), colnames(refGraph)))
        stop("Row and column names of 'refGraph' should be the same\n")
    }
  }

  if (!identical(rownames(measurementsMatrix), rownames(refGraph)))
    stop("Row and column names in 'measurementsMatrix' and 'refGraph' should be the same\n")

  if (fcOutput.na.rm && is.null(fcOutput))
    stop("When 'fcOutput.na.rm=TRUE then 'fcOutput' should be set.\n")

  var.names <- rownames(measurementsMatrix)
  n.var <- nrow(measurementsMatrix)

  ## by now we have to coerce it to a regular matrix
  ## but hopefully in the near future we can do [<- as well on a dspMatrix
  measurementsMatrix <- as.matrix(measurementsMatrix)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("'pairup.i' and 'pairup.j' should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (any(is.na(match(c(pairup.i, pairup.j), var.names))))
      warning("Some of the variables in 'pairup.i' or 'pairup.j' do not form part of 'measurementsMatrix' or 'refGraph'\n")

    pairup.i <- match(pairup.i, var.names)
    pairup.i <- pairup.i[!is.na(pairup.i)]

    pairup.j <- match(pairup.j, var.names)
    pairup.j <- pairup.j[!is.na(pairup.j)]

    pairup.ij.int <- intersect(pairup.i, pairup.j)
    pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
    pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

    nomeasurementsMask <- matrix(FALSE, nrow=n.var, ncol=n.var)
    nomeasurementsMask[as.matrix(
                       expand.grid(pairup.ij.int,
                                   c(pairup.i.noint, pairup.j.noint)))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.i.noint, pairup.j.noint))] <- TRUE
    nomeasurementsMask[as.matrix(expand.grid(pairup.ij.int, pairup.ij.int))] <- TRUE
    diag(nomeasurementsMask) <- FALSE
    nomeasurementsMask <- nomeasurementsMask | t(nomeasurementsMask)
    nomeasurementsMask <- !nomeasurementsMask
    measurementsMatrix[nomeasurementsMask] <- NA
  }

  upperTriRow <- row(measurementsMatrix)[upper.tri(measurementsMatrix) &
                                         !is.na(measurementsMatrix) & refGraph]
  upperTriCol <- col(measurementsMatrix)[upper.tri(measurementsMatrix) &
                                         !is.na(measurementsMatrix) & refGraph]
  measurementsUpperTriMatrix <- measurementsMatrix[upper.tri(measurementsMatrix) &
                                                   !is.na(measurementsMatrix) & refGraph]
  orderedMeasurements <- sort(measurementsUpperTriMatrix,index.return=TRUE,
                              decreasing=decreasing)

  if (!fcOutput.na.rm) {
    if (n == Inf)
      n <- length(measurementsUpperTriMatrix)

    edgeRnk <- data.frame(i=var.names[upperTriRow[orderedMeasurements$ix[1:n]]],
                          j=var.names[upperTriCol[orderedMeasurements$ix[1:n]]],
                          x=orderedMeasurements$x[1:n], stringsAsFactors=FALSE)
  } else {
    n2 <- length(measurementsUpperTriMatrix)

    edgeRnk <- data.frame(i=var.names[upperTriRow[orderedMeasurements$ix[1:n2]]],
                          j=var.names[upperTriCol[orderedMeasurements$ix[1:n2]]],
                          x=orderedMeasurements$x[1:n2], stringsAsFactors=FALSE)
  }

  if (!is.null(pairup.i)) {
    swapMask <- is.na(match(edgeRnk$i, var.names[pairup.i]))
    subRnk <- edgeRnk[swapMask, ]
    subRnk <- data.frame(subRnk$j, subRnk$i, subRnk$x, stringsAsFactors=FALSE)
    edgeRnk[swapMask, ] <- subRnk
  }

  if (!is.null(annotation)) {
    haveSYMBOL <- !is.na(match(paste0(gsub(".db", "", annotation), "SYMBOL"),
                               ls(sprintf("package:%s", annotation))))
    haveGNAME <- !is.na(match(paste0(gsub(".db", "", annotation), "GENENAME"),
                              ls(sprintf("package:%s", annotation))))
    if (!haveSYMBOL && !haveGNAME)
      stop("No SYMBOL nor GENENAME mappings for the annotation package %s\n", annotation)

    syms <- NA
    if (haveSYMBOL)
      syms <- unlist(AnnotationDbi::mget(unique(c(edgeRnk$i, edgeRnk$j)),
                       annotate::getAnnMap(map="SYMBOL", chip=annotation, type="db"),
                       ifnotfound=NA))
    else
      syms <- unlist(AnnotationDbi::mget(unique(c(edgeRnk$i, edgeRnk$j)),
                       annotate::getAnnMap(map="GENENAME", chip=annotation, type="db"),
                       ifnotfound=NA))
    edgeRnk$iSymbol=as.character(as.vector(syms[edgeRnk$i]))
    edgeRnk$jSymbol=as.character(as.vector(syms[edgeRnk$j]))
    edgeRnk <- edgeRnk[, c(1,2,4,5,3)]
    edgeRnk$iSymbol[is.na(edgeRnk$iSymbol)] <- as.vector(edgeRnk$i[is.na(edgeRnk$iSymbol)])
    edgeRnk$jSymbol[is.na(edgeRnk$jSymbol)] <- as.vector(edgeRnk$j[is.na(edgeRnk$jSymbol)])
  }

  if (!is.null(digits))
    edgeRnk$x <- round(edgeRnk$x, digits=digits)

  if (!haveMeasurements)
    edgeRnk <- edgeRnk[, -dim(edgeRnk)[2]]

  if (!is.null(fcOutput)) {
    if (class(fcOutput) != "list")
      stop("'fcOutput' should be the output of 'qpFunctionalCoherence'.\n")

    if (!all(names(fcOutput) == c("txRegNet", "txRegNetGO", "functionalCoherenceValues")))
      stop("'fcOutput' should be the output of 'qpFunctionalCoherence'.\n")

    if (is.null(pairup.i))
      stop("When 'fcOutput' is set 'pairup.i' and 'pairup.j' should be also set.\n")

    edgeRnk <- cbind(edgeRnk, funCoherence=NA_real_)
    fcMask <- !is.na(match(edgeRnk$i, names(fcOutput$functionalCoherenceValues)))
    if (any(fcMask)) {
      fcMask[fcMask] <- apply(edgeRnk[fcMask, ], 1, function(x, tx) !is.na(match(x[2], tx[[x[1]]])),
                              fcOutput$txRegNet)
      edgeRnk$funCoherence[fcMask] <- fcOutput$functionalCoherenceValues[edgeRnk$i[fcMask]]
    }

    if (!is.null(digits))
      edgeRnk$funCoherence <- round(edgeRnk$funCoherence, digits=digits)

    if (fcOutput.na.rm) {
      if (!any(!is.na(edgeRnk$funCoherence)))
        stop("No functional coherence value is different from NA. Please set 'fcOutput.na.rm=FALSE'\n")

      edgeRnk <- edgeRnk[!is.na(edgeRnk$funCoherence), ]

      if (n == Inf)
        n <- length(measurementsUpperTriMatrix)

      edgeRnk <- edgeRnk[1:n, ]
    }
  }

  rownames(edgeRnk) <- 1:dim(edgeRnk)[1]

  if (is.null(file))
    return(edgeRnk)
  else
    write.table(edgeRnk, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

  invisible(edgeRnk)
}



## function: qpPlotNetwork
## purpose: plot a given network using Rgraphviz
## parameters: g - the network to plot given as a graph
##             vertexSubset - a subset of vertices inducing the subgraph we want to plot
##             boundary - flag set to TRUE if we want also to plot the vertices connected
##                        to the given vertex subset including their connecting edges;
##                        FALSE (default) otherwise
##             minimumSizeConnComp - minimum size of the connected components we want to plot
##             pairup.i - subset of vertices to pair up with subset pairup.j
##             pairup.j - subset of vertices to pair up with subset pairup.i
##             highlight - subset of vertices to highlight
##             annotation - name of an annotation package to transform gene identifiers
##                          into gene symbols
##             layout - layout employed by Rgraphviz

qpPlotNetwork <- function(g, vertexSubset=graph::nodes(g), boundary=FALSE,
                          minimumSizeConnComp=2, pairup.i=NULL, pairup.j=NULL,
                          highlight=NULL, annotation=NULL, layout=c("twopi", "dot", "neato", "circo", "fdp")) {
  ## require(graph)

  layout <- match.arg(layout)

  if (any(is.na(match(graph::nodes(g), vertexSubset)))) {
    vertexSubsetNoMatch <- vertexSubset[is.na(match(vertexSubset, graph::nodes(g)))]
    if (length(vertexSubsetNoMatch) > 0 && is.null(annotation))
      stop("Vertex names in 'vertexSubset' ", vertexSubsetNoMatch, " do not form part of the vertices in 'g' and 'annotation' is set to NULL.\n")

    if (length(vertexSubsetNoMatch) > 0 && !is.null(annotation)) {
      vertexSubsetNoMatchIDs <- unlist(AnnotationDbi::mget(vertexSubsetNoMatch,
                                         AnnotationDbi::revmap(annotate::getAnnMap(map="SYMBOL", chip=annotation, type="db")),
                                         ifnotfound=NA))
      if (any(is.na(vertexSubsetNoMatchIDs)))
        stop("Vertex names in 'vertexSubset' ", vertexSubsetNoMatch[is.na(vertexSubsetNoMatchIDs)], " do not form part of the vertices in 'g' and identifiers could not be found through the SYMBOL map from 'annotation'.\n")

      vertexSubset <- c(setdiff(vertexSubset, vertexSubsetNoMatch), vertexSubsetNoMatchIDs)
    }

    if (boundary) {
      bd <- boundary(vertexSubset, g)
      bd <- bd[sapply(bd, length) > 0]
      vertexSubset <- unique(c(vertexSubset, unlist(bd)))
    }
    g <- subGraph(vertexSubset, g)
    if (numEdges(g) == 0)
      stop("The resulting graph has no edges\n")
  }

  if (minimumSizeConnComp > 1) {
    gCc <- connComp(g)
    gCcOfMinSize <- gCc[sapply(gCc, length) >= minimumSizeConnComp]
    vertexSubset <- unique(unlist(gCcOfMinSize))

    if (any(is.na(match(graph::nodes(g), vertexSubset)))) {
      g <- subGraph(vertexSubset, g)
      if (numEdges(g) == 0)
        stop("The resulting graph has no edges\n")
    }
  }

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("'pairup.i' and 'pairup.j' should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (sum(is.na(match(graph::nodes(g), c(pairup.i, pairup.j)))) > 0)
      warning("Some of the vertices in the resulting graph do not form part of 'pairup.i' nor 'pairup.j'\n")

    edL <- matrix(unlist(sapply(graph::nodes(g),
                                function(x) t(cbind(x, graph::edges(g)[[x]])), USE.NAMES=FALSE)),
                   ncol=2, byrow=TRUE)
    edL <- unique(t(apply(edL, 1, sort)))
    mask <- apply(edL, 1, function(x) sum(!is.na(match(x, pairup.i)))*sum(!is.na(match(x, pairup.j)))) == 0
    if (sum(mask) > 0)
      g <- removeEdge(from=edL[mask, 1], to=edL[mask, 2], g)

    edgemode(g) <- "directed"
    g.iNodes <- graph::nodes(g)[!is.na(match(graph::nodes(g), pairup.i))]
    wrongEdges <- boundary(setdiff(graph::nodes(g), g.iNodes), g)
    wrongEdges <- wrongEdges[sapply(wrongEdges, length) > 0]
    wrongEdges <- matrix(unlist(sapply(names(wrongEdges), function(x) t(cbind(x, wrongEdges[[x]])), USE.NAMES=FALSE)),
                         ncol=2, byrow=TRUE)
    g <- removeEdge(from=wrongEdges[, 1], to=wrongEdges[, 2], g)

    nodeLabels <- graph::nodes(g)
    if (!is.null(annotation)) {
      haveSYMBOL <- !is.na(match(paste0(gsub(".db", "", annotation), "SYMBOL"),
                                 ls(sprintf("package:%s", annotation))))
      haveGNAME <- !is.na(match(paste0(gsub(".db", "", annotation), "GENENAME"),
                                ls(sprintf("package:%s", annotation))))
      if (!haveSYMBOL && !haveGNAME)
        stop("No SYMBOL nor GENENAME mappings for the annotation package %s\n", annotation)

      if (haveSYMBOL)
        nodeLabels <- unlist(AnnotationDbi::mget(graph::nodes(g),
                               annotate::getAnnMap(map="SYMBOL", chip=annotation, type="db"),
                               ifnotfound=NA))
      else
        nodeLabels <- unlist(AnnotationDbi::mget(graph::nodes(g),
                               annotate::getAnnMap(map="GENENAME", chip=annotation, type="db"),
                               ifnotfound=NA))
      nodeLabels[is.na(nodeLabels)] <- graph::nodes(g)[is.na(nodeLabels)]
    }

    pkg <- "Rgraphviz"
    if (!library(pkg, logical.return=TRUE, character.only=TRUE)) {
      warning("qpPlotNetwork() requires package 'Rgraphviz' to plot the network and does not seem to be installed\n")
      return(invisible(g))
    }

    mkNodeAttrs <- get("makeNodeAttrs", mode="function")
    nodattr <- mkNodeAttrs(g, label=nodeLabels, shape="ellipse", fixedsize=FALSE,
                           fillcolor=grey(0.9), fontcolor="black")
    nodattr$fillcolor[g.iNodes] <- grey(0.65)
  } else {
    nodeLabels <- graph::nodes(g)
    if (!is.null(annotation)) {
      haveSYMBOL <- !is.na(match(paste0(gsub(".db", "", annotation), "SYMBOL"),
                                 ls(sprintf("package:%s", annotation))))
      haveGNAME <- !is.na(match(paste0(gsub(".db", "", annotation), "GENENAME"),
                                ls(sprintf("package:%s", annotation))))
      if (!haveSYMBOL && !haveGNAME)
        stop("No SYMBOL nor GENENAME mappings for the annotation package %s\n", annotation)

      if (haveSYMBOL)
        nodeLabels <- unlist(AnnotationDbi::mget(graph::nodes(g),
                               annotate::getAnnMap(map="SYMBOL", chip=annotation, type="db"),
                               ifnotfound=NA))
      else
        nodeLabels <- unlist(AnnotationDbi::mget(graph::nodes(g),
                               annotate::getAnnMap(map="GENENAME", chip=annotation, type="db"),
                               ifnotfound=NA))
      nodeLabels[is.na(nodeLabels)] <- graph::nodes(g)[is.na(nodeLabels)]
    }

    pkg <- "Rgraphviz"
    if (!library(pkg, logical.return=TRUE, character.only=TRUE)) {
      warning("qpPlotNetwork() requires package 'Rgraphviz' to plot the network and does not seem to be installed\n")
      return(invisible(g))
    }

    mkNodeAttrs <- get("makeNodeAttrs", mode="function")
    nodattr <- mkNodeAttrs(g, label=nodeLabels, shape="ellipse", fixedsize=FALSE,
                           fillcolor=grey(0.9), fontcolor="black")
  }

  if (!is.null(highlight)) {
    g.hNodes <- graph::nodes(g)[!is.na(match(graph::nodes(g), highlight))]
    nodattr$fontcolor[g.hNodes] <- "red"
  }

  plot(g, layout, nodeAttrs=nodattr, lwd=2)

  invisible(g)
}



## function: qpPlotMap
## purpose: plot a map of associated pairs defined by adjusted p-values
## parameters: p.valueMatrix - a squared symmetric matrix with raw p-values for all pairs
##             markerPos - a two-column matrix containing chromosome and position of each feature
##             genePos - a two-column matrix containing chromosome and position of each feature
##             chrLen - a named vector with chromosome lengths
##             p.value - adjusted p-value cutoff
##             adjust.method - method employed to adjust the raw p-values. It is passed in a call
##                             call to p.adjust() in its 'method' argument
##             xlab - label for the x-axis
##             ylab - label for the y-axis
##             main - main title of the plot, set to the empty string by default

qpPlotMap <- function(p.valueMatrix, markerPos, genePos, chrLen,
                      p.value=0.05, adjust.method="holm",
                      xlab="Ordered Markers", ylab="Ordered Genes",
                      main="", ...) {

  ## get all variables
  var.names <- rownames(p.valueMatrix)

  ## get positions ordered by chr
  markerPos <- markerPos[order(markerPos[, 1], markerPos[,2]), ]
  genePos <- genePos[order(genePos[, 1], genePos[,2]), ]

  ## re-scale genes and marker positions between 0 and 1
  chrRelLen <- chrLen/sum(chrLen)
  chrRelCumLen <- c(0, cumsum(chrRelLen)[-length(chrRelLen)])
  geneRelPos <- chrRelCumLen[genePos[ ,1]] +
                chrRelLen[genePos[ ,1]]*(genePos[ ,2]/chrLen[genePos[ ,1]]) 
  names(geneRelPos) <- rownames(genePos)
  markerRelPos <- chrRelCumLen[markerPos[ ,1]] +
                  chrRelLen[markerPos[ ,1]]*(markerPos[ ,2]/chrLen[markerPos[ ,1]])
  names(markerRelPos) <- rownames(markerPos)

  ## build map matrix
  tp <- qpTopPairs(p.valueMatrix,
                   pairup.i=intersect(names(markerRelPos), var.names),
                   pairup.j=intersect(names(geneRelPos), var.names), n=Inf)
  tp <- cbind(tp, padj=p.adjust(tp$x, method=adjust.method))
  tp <- tp[tp$padj < p.value, ]
  if (nrow(tp) == 0)
    stop(sprintf("No association meets the p.value cutoff at %.2f. Increase the p.value or change the p-value adjusting method\n", p.value))

  tp <- tp[(tp$i %in% names(markerRelPos)) & (tp$j %in% names(geneRelPos)), ]
  colnames(tp) <- c("i", "j", "pval", "adjpval")

  ## jitter a bit equal positions
  lastIdx <- length(markerRelPos)
  while(sum(markerRelPos - c(0, markerRelPos[-lastIdx]) == 0) > 0) {
    wh.equal <- which(markerRelPos - c(0, markerRelPos[-lastIdx]) == 0)
    markerRelPos[wh.equal] <- markerRelPos[wh.equal] +
                            (markerRelPos[wh.equal + 1] - markerRelPos[wh.equal]) / 2
  }

  lastIdx <- length(geneRelPos)
  while(sum(geneRelPos - c(0, geneRelPos[-lastIdx]) == 0) > 0) {
    wh.equal <- which(geneRelPos - c(0, geneRelPos[-lastIdx]) == 0)
    geneRelPos[wh.equal] <- geneRelPos[wh.equal] +
                            (geneRelPos[wh.equal + 1] - geneRelPos[wh.equal]) / 2
  }

  ## plot the map
  x <- markerRelPos[tp[, 1]]
  y <- geneRelPos[tp[, 2]]
  plot(x, y, pch=".", xlab=xlab, ylab=ylab, main=main,
       axes=FALSE, xlim=range(x), ylim=range(y), ...)
  segments(chrRelCumLen, 0, chrRelCumLen, 1, col=gray(0.75), lty="dotted", lwd=2)
  segments(0, chrRelCumLen, 1, chrRelCumLen, col=gray(0.75), lty="dotted", lwd=2)
  axis(1, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
       labels=names(chrLen), tick=FALSE, cex.axis=0.9)
  axis(2, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
       labels=names(chrLen), tick=FALSE, cex.axis=0.9, las=1)

  invisible(tp)
}

#####################
# PRIVATE FUNCTIONS #
#####################

## function: qpEdgePACSE
## purpose: calculate the standard errors for the partial correlations of the
##          edges of an undirected graphical Gaussian Markov model according to
##          the method by:
##
##          Roverato and Whittaker. Standard errors for the parameters of
##          graphical Gaussian models, STATISTICS AND COMPUTING, 6:297-302, 1996)

## parameters: S - estimate of the sample covariance matrix
##             A - adjacency matrix of the graph and thus it is assumed that the diagonal
##                 is set to either 0s or FALSE truth values since there should be no loops
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two members: K - the concentration matrix; SE the matrix
##         with the standard errors of the edges

.qpEdgePACSE <- function(S, A, R.code.only=FALSE) {

  if (!R.code.only) {
    ## this should change so that the entire algorithm deals with *Matrix classes from the Matrix package
    return(qpgraph:::.qpFastPACSE(as(S, "matrix"), as(A, "matrix")));
  }

  A <- as(A, "matrix") ## idem

  n.var <- nrow(A)

  A <- A + diag(n.var) ## in the code below we need 1s in the main diagonal and
                       ## then at the same time we make sure we get a 0-1 matrix
                       ## as a truth value + 0 or 1 equals a number

  A[col(A) > row(A)] <- NA

  # selection row and column indices corr. to the non-zero elem.
  A[A == 0] <- NA
  r.nz <- c(row(A))[!is.na(A)]
  c.nz <- c(col(A))[!is.na(A)]

  # computation of the Fisher information matrix
  Iss <- S[c.nz,c.nz] * S[r.nz,r.nz] +
         S[c.nz,r.nz] * S[r.nz,c.nz]
  Iss <- solve(Iss)
  IssI <- matrix(rep(0, length(r.nz) * length(c.nz)), nrow=length(r.nz))
  diag(IssI) <- 1
  IssI[cbind((1:length(r.nz))[r.nz==c.nz], (1:length(r.nz))[r.nz==c.nz])] <- 2

  FISHER <- IssI %*% Iss %*% IssI

  # standard errors are in the diagonal of the Fisher information matrix
  FSHR <- diag(FISHER)
  SE <- matrix(NA, nrow(A), nrow(A))
  SE[cbind(r.nz,c.nz)] <- SE[cbind(c.nz,r.nz)] <- FSHR
  diag(SE) <- NA

  return(SE)
}



## function: qpIPFpass
## purpose: implement the Iterative Proportional Fitting (IPF) algorithm. Part of
##          this R code has been borrowed from an implementation by Graham Wills
##          in June of 1992
## parameters: Vf - matrix to adjust
##             Vn - matrix to adjust
##             clqlst - list of (maximal) cliques
## return: the input matrix adjusted to the constraints of the list of cliques

.qpIPFpass <- function(Vf, Vn, clqlst) {
  n.var <- nrow(Vf)
  firstclq <- 1

  if (clqlst[[1]][1] > n.var) { # if the clique list has vertex-clique indices
    firstclq <- n.var + 1       # at the beginning
  }

  V <- Vn
  for(i in firstclq:length(clqlst)) {
    V <- qpgraph:::.qpIPFstep(Vf, V, i, clqlst)
  }

  return(V)
}



## function: qpIPFstep
## purpose: implement the Iterative Proportional Fitting (IPF) algorithm. Part of
##          this R code has been borrowed from an implementation by Graham Wills
##          in June of 1992
## parameters: Vf - matrix to adjust
##             Vn - matrix to adjust
##             wh - clique index
##             clqlst - (maximal) clique
## return: the input matrix adjusted to the constraints of the clique

.qpIPFstep <- function(Vf, Vn, wh, clqlst) {
  a <- clqlst[[wh]]
  b <- (1:length(Vf[, 1]))[ - a]
  Vfaa <- Vf[a, a]
  Vni <- solve(Vn[a, a])
  Bnba <- Vn[b, a] %*% Vni
  Vnbba <- Vn[b, b] - Vn[b, a] %*% Vni %*% Vn[a, b]
  V <- Vf
  V[b, a] <- Bnba %*% Vfaa
  V[a, b] <- t(V[b, a])
  V[b, b] <- Vnbba + Bnba %*% Vfaa %*% t(Bnba)

  return(V)
}



## function: qpFilterByGO
## purpose: filter an input vector of Entrez Gene IDs returning only those that
##          have Gene Ontology (GO) annotations on a specified ontology branch
## parameters: entrezGeneIds - Entrez Gene Identifiers to be filtered
##             chip - .db package name containing the GO annotations
##             ontologyType - either "BP", or "MF" or "CC"
## return: the subset entrezGeneIds for which GO annotations of the specified
##         ontology branch are found

.qpFilterByGO <- function(entrezGeneIds, chip, ontologyType=c("BP", "MF", "CC")) {
  ontologyType <- match.arg(ontologyType)

  if (length(entrezGeneIds) == 0)
    stop("qpgraph:::.qpFilterByGO: no input Entrez gene identifiers\n")

  haveGo <- sapply(AnnotationDbi::mget(entrezGeneIds,
                                       annotate::getAnnMap(map="GO", chip=chip, type="db"),
                                       ifnotfound=NA),
                   function(x) {
                     if (length(x) == 1 && is.na(x))
                       FALSE
                     else {
                       onts <- subListExtract(x, "Ontology", simplify=TRUE)
                       ontologyType %in% onts
                     }
                   })

  filteredIds <- names(haveGo)[haveGo]

  return(filteredIds)
}


## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: qpIsPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.qpIsPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
 (paste("package:", name, sep="") %in% search()) ||
 (name %in% loadedNamespaces())
}    


## function: qpIsPackageInstalled
## purpose: to check whether the package specified by the name given in
##          the input argument is installed.
## parameters: name - package name
## return: TRUE if the package is installed, FALSE otherwise

.qpIsPackageInstalled <- function(name) {
  ## Purpose: is package 'name' installed?
  ## --------------------------------------------------
 (name %in% rownames(installed.packages()))
}


## function: clPrCall
## purpose: is a copy of the function clusterCall() from the 'snow' package
##          but allows the slave loops to report progress which is then
##          reported to the console of the master node
## parameters: cl - cluster data from makeCluster()
##             fun - function to call at each node
##             n.adj - total number of adjacencies through which computations
##                     should be made
##             ... - parameters passed to the function specified at 'fun'
## return: the result just as clusterCall() would do

clPrCall <- function(cl, fun, n.adj, ...) {
  ## checkCl <- get("checkCluster", mode="function")
  ## sndCall <- get("sendCall", mode="function")
  ## rcv1Result <- get("recvOneResult", mode="function")
  ## check4RmtErrors <- get("checkForRemoteErrors", mode="function")

  snow::checkCluster(cl)
  for (i in seq(along = cl))
    snow::sendCall(cl[[i]], fun, list(...))

  i <- rep(0, length(cl))
  k <- 0
  ppct <- -1
  pb <- txtProgressBar(style=3)

  r <- vector(mode="list", length=length(cl))
  foundError <- FALSE
  nodesDone <- 0
  while (nodesDone < length(cl) && !foundError) {
    r1 <- snow::recvOneResult(cl)
    if (!is.null(r1$tag)) {
      ## message("received value ", r1$value, " from node ", r1$node, " with tag ", r1$tag)
      if (r1$tag != "UPDATE") {
        if (inherits(r1, "try-error")) {
          stop("at least one node produced an error: ", r)
          foundError <- TRUE
        }
      } else {
        k <- k - i[r1$node] + r1$value
        i[r1$node] <- r1$value
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
    } else {
      r1 <- snow::checkForRemoteErrors(r1)
      r[[r1$node]] <- r1$value
      nodesDone <- nodesDone + 1
    }
  }

  close(pb)

  r
}

##########################################################################
## PRIVATE FUNCTIONS THAT ARE ENTRY POINTS TO THE C CODE OF THE PACKAGE ##
##########################################################################

.qpFastNrr <- function(X, I, Y, q, restrict.Q, fix.Q, nTests, alpha,
                       pairup.i.noint, pairup.j.noint, pairup.ij.int,
                       exact.test, verbose, startTime, nAdj2estimateTime) {

  ## nLevels <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
  nLevels <- rep(NA_integer_, times=ncol(X))
  nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
  if (any(nLevels[I] == 1))
    stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

  return(.Call("qp_fast_nrr", X, as.integer(I), as.integer(nLevels),
                             as.integer(Y), as.integer(q), restrict.Q, ## restrict.Q can be a matrix
                             as.integer(fix.Q), as.integer(nTests), as.double(alpha),
                             as.integer(pairup.i.noint), as.integer(pairup.j.noint),
                             as.integer(pairup.ij.int), as.integer(exact.test),
                             as.integer(verbose), as.double(startTime),
                             as.integer(nAdj2estimateTime), .GlobalEnv))
}

.qpFastNrrIdenticalQs <- function(X, q, restrict.Q, fix.Q, nTests, alpha,
                                  pairup.i.noint, pairup.j.noint, pairup.ij.int,
                                  verbose, startTime, nAdj2estimateTime) {

  return(.Call("qp_fast_nrr_identicalQs", X, as.integer(q), as.integer(restrict.Q),
                                         as.integer(fix.Q),
                                         as.integer(nTests), as.double(alpha),
                                         as.integer(pairup.i.noint),
                                         as.integer(pairup.j.noint),
                                         as.integer(pairup.ij.int),
                                         as.integer(verbose), as.double(startTime),
                                         as.integer(nAdj2estimateTime), .GlobalEnv))
}

.qpFastNrrPar <- function(X, I, Y, q, restrict.Q, fix.Q, nTests, alpha,
                          pairup.i.noint, pairup.j.noint, pairup.ij.int,
                          exact.test, verbose, estimateTime, nAdj2estimateTime) {
  clOpt <- get("getClusterOption", mode="function")
  myMaster <- clOpt("masterNode")

  startTime <- 0
  if (estimateTime)
    startTime <- proc.time()["elapsed"]

  ## nLevels <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
  nLevels <- rep(NA_integer_, times=ncol(X))
  nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
  if (any(nLevels[I] == 1))
    stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

  ## clusterRank and clusterSize should have been defined by the master node
  return(.Call("qp_fast_nrr_par", X, as.integer(I), as.integer(nLevels),
                                 as.integer(Y), as.integer(q), restrict.Q, ## restrict.Q can be a matrix
                                 as.integer(fix.Q), as.integer(nTests), as.double(alpha),
                                 as.integer(pairup.i.noint), as.integer(pairup.j.noint),
                                 as.integer(pairup.ij.int), as.integer(exact.test),
                                 as.integer(verbose), as.double(startTime),
                                 as.integer(nAdj2estimateTime), as.integer(get("clusterRank")),
                                 as.integer(get("clusterSize")), myMaster, .GlobalEnv))
}

.qpFastNrrIdenticalQsPar <- function(X, q, restrict.Q, fix.Q, nTests, alpha,
                                     pairup.i.noint, pairup.j.noint,
                                     pairup.ij.int, verbose, estimateTime,
                                     nAdj2estimateTime) {
  clOpt <- get("getClusterOption", mode="function")
  myMaster <- clOpt("masterNode")

  startTime <- 0
  if (estimateTime)
    startTime <- proc.time()["elapsed"]

  ## clusterRank and clusterSize should have been defined by the master node
  return(.Call("qp_fast_nrr_identicalQs_par",X,as.integer(q), as.integer(restrict.Q),
                                             as.integer(fix.Q), as.integer(nTests),
                                             as.double(alpha), as.integer(pairup.i.noint),
                                             as.integer(pairup.j.noint),
                                             as.integer(pairup.ij.int),
                                             as.integer(verbose), as.double(startTime),
                                             as.integer(nAdj2estimateTime),
                                             as.integer(get("clusterRank")),
                                             as.integer(get("clusterSize")),
                                             myMaster, .GlobalEnv))
}

.qpFastEdgeNrr <- function(X, S, i, j, q, restrict.Q, fix.Q, nTests, alpha) {
  Sx <- n <- p <- NULL
  if (!is.null(S)) {
    Sx <- S@ssd@x
    n <- S@n
    p <- ncol(S)
  } else {
    p <- ncol(X)
    n <- nrow(X)
  }

  return(.Call("qp_fast_edge_nrr", X, Sx, p, n, as.integer(i), as.integer(j),
                                  as.integer(q), as.integer(restrict.Q), as.integer(fix.Q),
                                  as.integer(nTests), as.double(alpha)))
}

.qpFastEdgeNrrHMGM <- function(X, I, nLevels, Y, ssd, mapX2ssd, i, j, q, restrict.Q,
                               fix.Q, nTests, alpha, exact.test) {
  ssdx <- NULL
  if (!is.null(ssd)) {
    if (class(ssd) != "SsdMatrix")
      stop("qpgraph:::.qpFastEdgeNrrHMGM: the ssd argument should be an object of class SsdMatrix\n")
    ssdx <- ssd@ssd@x
  }

  return(.Call("qp_fast_edge_nrr_hmgm", X, as.integer(I), as.integer(nLevels),
                                        as.integer(Y), ssdx, as.integer(mapX2ssd),
                                        as.integer(i),as.integer(j), as.integer(q),
                                        as.integer(restrict.Q), as.integer(fix.Q),
                                        as.integer(nTests), as.double(alpha), as.integer(exact.test)))
}

.qpFastEdgeNrrHMGMsml <- function(X, cumsum_sByChr, s, gLevels, XEP, I, nLevels, Y,
                                  ssd, mapX2ssd, i, j, q, restrict.Q, fix.Q,
                                  nTests, alpha, exact.test) {

  ssdx <- NULL
  if (!is.null(ssd)) {
    if (class(ssd) != "SsdMatrix")
      stop("qpgraph:::.qpFastEdgeNrrHMGM: the ssd argument should be an object of class SsdMatrix\n")
    ssdx <- ssd@ssd@x
  }

  return(.Call("qp_fast_edge_nrr_hmgm_sml", X, as.integer(cumsum_sByChr), as.integer(s),
                                            as.integer(gLevels), XEP, as.integer(I),
                                            as.integer(nLevels), as.integer(Y), ssdx,
                                            as.integer(mapX2ssd), as.integer(i),
                                            as.integer(j), as.integer(q),
                                            as.integer(restrict.Q), as.integer(fix.Q),
                                            as.integer(nTests), as.double(alpha),
                                            as.integer(exact.test)))
}

.qpFastCliquerGetCliques <- function(A, clqspervtx, verbose) {
  return(.Call("qp_fast_cliquer_get_cliques", A, clqspervtx, verbose))
}

.qpFastUpdateCliquesRemoving <- function(A, clqlst, v, w, verbose) {
  return(.Call("qp_fast_update_cliques_removing", A, clqlst, v, w, verbose))
}

.qpFastPACSE <- function(Shat, A) {
  return(.Call("qp_fast_pac_se", Shat, A))
}

.qpFastIPF <- function(vv, clqlst, tol = 0.001, verbose = FALSE) {
  return(.Call("qp_fast_ipf", vv, clqlst, tol, verbose))
}

.qpFastHTF <- function(S, A, tol = 0.001, verbose = FALSE) {
  return(.Call("qp_fast_htf", S, A, tol, verbose))
}

.qpCliqueNumberLowerBound <- function(A, return.vertices, approx.iter, verbose) {
 return(.Call("qp_clique_number_lb", A, return.vertices, as.integer(approx.iter), verbose))
}

.qpCliqueNumberOstergard <- function(A, return.vertices, verbose) {
 return(.Call("qp_clique_number_os", A, return.vertices, verbose))
}

qpCov <- function(X, corrected=TRUE) {
  ## return(new("dspMatrix", Dim=c(ncol(X), ncol(X)),
  ##            Dimnames=list(colnames(X), colnames(X)),
  ##            x = .Call("qp_cov_upper_triangular", X, as.integer(corrected))))
  return(.Call("qp_fast_cov_upper_triangular", X, as.integer(corrected)))
}

