## qpgraph package - this R code implements functions to learn qp-graphs from
##                   data, to estimate Pearson and partial correlations and
##                   to interact with microarray data in order to build network
##                   models of molecular regulation
##
## Copyright (C) 2008 R. Castelo and A. Roverato
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
## parameters: data - data set from where to estimate the non-rejection rates
##             q - partial-correlation order to be employed
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
##             R.code.only - flag set to FALSE when using the C implementation
## return: a matrix with the estimates of the non-rejection rates

setGeneric("qpNrr", function(data, ...) standardGeneric("qpNrr"))

# data comes as an ExpressionSet object
setMethod("qpNrr", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpNrr(S, N, ...)
          })

# data comes as a data frame
setMethod("qpNrr", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpNrr(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpNrr", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            S <- cov(data)
            N <- length(data[,1])
            if (!is.null(colnames(data))) 
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpNrr(S, N, ...)
          })

.qpNrr <- function(S, N, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                   pairup.j=NULL, verbose=TRUE, R.code.only=FALSE) {

  var.names <- rownames(S)
  n.var <- nrow(S)

  if (alpha < 0.0 || alpha > 1.0) {
    stop(sprintf("significance level alpha is %.2f and it must lie in the interval [0,1]\n",alpha))
  }

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j must both either be set to NULL or contain subsets of variables\n")

  if (is.null(pairup.i))
    pairup.i <- 1:n.var
  else {
    pairup.i <- match(pairup.i, var.names)
    if (sum(is.na(pairup.i)) > 0)
      stop("pairup.i is not a subset of the variables forming the data\n")
  }

  if (is.null(pairup.j))
    pairup.j <- 1:n.var
  else {
    pairup.j <- match(pairup.j, var.names)
    if (sum(is.na(pairup.j)) > 0)
      stop("pairup.j is not a subset of the variables forming the data\n")
  }

  # pair the two sets pairup.i and pairup.j without pairing the same variable
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

  if (!R.code.only) {
    nrrMatrix <- qpgraph:::.qpFastNrr(S, N, q, nTests, alpha, pairup.i.noint,
                                 pairup.j.noint, pairup.ij.int, verbose)
    rownames(nrrMatrix) <- colnames(nrrMatrix) <- var.names

    return(nrrMatrix)
  }

  if (q > n.var - 2)
    stop(paste("q=",q," > n.var-2=",n.var-2))

  if (q < 0)
    stop(paste("q=",q," < 0"))

  if (q > N - 3)
    stop(paste("q=",q," > N-3=",N-3))

  nrrMatrix <- matrix(as.double(NA), n.var, n.var)
  rownames(nrrMatrix) <- colnames(nrrMatrix) <- var.names
  ppct <- -1
  k <- 0

  # intersection variables against ij-exclusive variables
  for (i in pairup.ij.int) {
    for (j in c(pairup.i.noint,pairup.j.noint)) {
      nrrMatrix[j,i] <- nrrMatrix[i,j] <-
        qpEdgeNrr(S, N, i, j, q, nTests, alpha, R.code.only=TRUE)
      k <- k + 1
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose) {
        if (pct %% 10 == 0) {
          cat(pct)
        } else {
          cat(".")
        }
        ppct <- pct
      }
    }
  }

  # i-exclusive variables against j-exclusive variables
  for (i in pairup.i.noint) {
    for (j in pairup.j.noint) {
      nrrMatrix[j,i] <- nrrMatrix[i,j] <-
        qpEdgeNrr(S, N, i, j, q, nTests, alpha, R.code.only=TRUE)
      k <- k + 1
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose) {
        if (pct %% 10 == 0) {
          cat(pct)
        } else {
          cat(".")
        }
        ppct <- pct
      }
    }
  }

  # intersection variables against themselves (avoiding pairing of the same)
  for (i in 1:(l.int-1)) {
    i2 <- pairup.ij.int[i]

    for (j in (i+1):l.int) {
      j2 <- pairup.ij.int[j]
      nrrMatrix[j2,i2] <- nrrMatrix[i2,j2] <-
        qpEdgeNrr(S, N, i2, j2, q, nTests, alpha, R.code.only=TRUE)
      k <- k + 1
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose) {
        if (pct %% 10 == 0) {
          cat(pct)
        } else {
          cat(".")
        }
        ppct <- pct
      }
    }
  }

  if (verbose) {
    cat("\n")
  }

  return(nrrMatrix)
}



## function: qpAvgNrr
## purpose: estimate average non-rejection rates for every pair of variables
## parameters: data - data set from where to estimate the average non-rejection
##                    rates
##             qOrders - either a number of partial-correlation orders or a
##                       vector of particular orders to be employed in the
##                       calculation
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
##             R.code.only - flag set to FALSE when using the C implementation
## return: a matrix with the estimates of the average non-rejection rates

setGeneric("qpAvgNrr", function(data, ...) standardGeneric("qpAvgNrr"))

# data comes as an ExpressionSet object
setMethod("qpAvgNrr", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpAvgNrr(S, N, ...)
          })

# data comes as a data frame
setMethod("qpAvgNrr", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpAvgNrr(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpAvgNrr", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            S <- cov(data)
            N <- length(data[,1])
            if (!is.null(colnames(data))) 
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpAvgNrr(S, N, ...)
          })

.qpAvgNrr <- function(S, N, qOrders=4, nTests=100, alpha=0.05, pairup.i=NULL,
                      pairup.j=NULL, type=c("arith.mean"), verbose=TRUE,
                      R.code.only=FALSE) {

  type <- match.arg(type)

  var.names <- rownames(S)
  n.var <- nrow(S)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j must both either be set to NULL or contain subsets of variables\n")

  if (length(qOrders) == 1) {
    if (qOrders > min(n.var, N) - 3)
      stop(sprintf("qOrders=%d is larger than the number of available q-orders for the given data set (%d)\n",
                   qOrders, min(n.var, N) - 3))

    qOrders <- as.integer(round(seq(1, min(n.var, N) - 3,
                                    by=(min(n.var, N) - 3) / qOrders), digits=0))
  } else {
    qOrders <- as.integer(qOrders)
    if (min(qOrders) < 1 || max(qOrders) > min(n.var,N))
      stop(sprintf("for the given data set q-orders must lie in the range [1,%d]\n",
                   min(n.var,N)))
  }

  w <- 1 / length(qOrders)
  avgNrrMatrix <- matrix(0,nrow=n.var,ncol=n.var)
  rownames(avgNrrMatrix) <- colnames(avgNrrMatrix) <- var.names
  for (q in qOrders) {
    if (verbose)
      cat(sprintf("q=%d\n",q))

    avgNrrMatrix <- avgNrrMatrix +
                    w * qpgraph:::.qpNrr(S, N, q, nTests, alpha, pairup.i,
                                         pairup.j, verbose, R.code.only)
  }

  return(avgNrrMatrix)
}



## function: qpEdgeNrr
## purpose: estimate the non-rejection rate for one pair of variables as the
##          fraction of tests that accept the null hypothesis of independence given
##          a set of randomly sampled q-order conditionals
## parameters: S - sample covariance matrix of the data
##             N - sample size
##             i - index in S (row or column) matching one of the two variables
##             j - index in S (row or column) matching the other variable
##             q - partial-correlation order
##             nTests - number of tests to perform
##             alpha - significance level of each test (Type-I error probability)
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             R.code.only - flag set to FALSE when using the C implementation
## return: an estimate of the non-rejection rate for the particular given pair of
##         variables

setGeneric("qpEdgeNrr", function(data, ...) standardGeneric("qpEdgeNrr"))

# data comes as an ExpressionSet object
setMethod("qpEdgeNrr", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpEdgeNrr(S, N, ...)
          })

# data comes as a data frame
setMethod("qpEdgeNrr", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpEdgeNrr(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpEdgeNrr", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            # if the matrix is squared let's assume then that it is the sample
            # covariance matrix and that the sample size is the next parameter
            if (nrow(data) != ncol(data)) {
              S <- cov(data)
              N <- length(data[,1])
              if (!is.null(colnames(data))) 
                rownames(S) <- colnames(S) <- 1:nrow(S)
              else
                rownames(S) <- colnames(S) <- colnames(data)
              qpgraph:::.qpEdgeNrr(S, N, ...)
            } else {
              S <- data
              qpgraph:::.qpEdgeNrr(S, ...)
            }
          })

.qpEdgeNrr <- function(S, N, i=1, j=2, q=1, nTests=100, alpha=0.05,
                       R.code.only=FALSE) {
  if (is.character(i)) {
    if (is.na(match(i, colnames(S))))
      stop(sprintf("i=%s does not form part of the variable names of the data\n",i))
    i <- match(i,colnames(S))
  }

  if (is.character(j)) {
    if (is.na(match(j, colnames(S))))
      stop(sprintf("j=%s does not form part of the variable names of the data\n",j))
    j <- match(j,colnames(S))
  }

  n.var  <- nrow(S)

  if (q > n.var-2)
    stop(paste("q=",q," > n.var-2=",n.var-2))

  if (q < 0)
    stop(paste("q=",q," < 0"))

  if (q > N-3)
    stop(paste("q=",q," > N-3=",N-3))

  if (!R.code.only) {
    return(qpgraph:::.qpFastEdgeNrr(S, N, i, j, q, nTests, alpha));
  }

  pop <- (1:n.var)[-c(i, j)]

  thr    <- qt(p=1-(alpha/2),df=N-q-2,lower.tail=TRUE,log.p=FALSE)
  lambda <- c()
  for (k in 1:nTests) {
    sp <- sample(pop, q, rep=F)
    cit <- qpgraph:::.qpCItest(S, N, i, j, sp, R.code.only=TRUE)
    lambda  <- c(lambda,abs(cit$t.value))
  }

  nAcceptedTests <- sum(lambda < thr)

  return(nAcceptedTests / nTests)
}



## function: qpCItest
## purpose: perform a conditional independence test between two variables given
##          a conditioning set
## parameters: data - data where to perform the test
##             N - sample size (when data is directly the sample covariance matrix)
##             i - index in S (row or column) matching one of the two variables
##             j - index in S (row or column) matching the other variable
##             Q - indexes in S of the variables forming the conditioning set
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two members, the t-statistic value and the p-value
##         on rejecting the null hypothesis of independence

setGeneric("qpCItest", function(data, ...) standardGeneric("qpCItest"))

# data comes as an ExpressionSet object
setMethod("qpCItest", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpCItest(S, N, ...)
          })

# data comes as a data frame
setMethod("qpCItest", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpCItest(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpCItest", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            # if the matrix is squared let's assume then that it is the sample
            # covariance matrix and that the sample size is the next parameter
            if (nrow(data) != ncol(data)) {
              S <- cov(data)
              N <- length(data[,1])
              if (!is.null(colnames(data))) 
                rownames(S) <- colnames(S) <- 1:nrow(S)
              else
                rownames(S) <- colnames(S) <- colnames(data)
              qpgraph:::.qpCItest(S, N, ...)
            } else {
              S <- data
              qpgraph:::.qpCItest(S, ...)
            }
          })

.qpCItest <- function(S, N, i=1, j=2, Q=c(), R.code.only=FALSE) {

  if (is.character(i)) {
    if (is.na(match(i, colnames(S))))
      stop(sprintf("i=%s does not form part of the variable names of the data\n",i))
    i <- match(i,colnames(S))
  }

  if (is.character(j)) {
    if (is.na(match(j, colnames(S))))
      stop(sprintf("j=%s does not form part of the variable names of the data\n",j))
    j <- match(j,colnames(S))
  }

  if (is.character(Q)) {
    if (sum(is.na(match(Q, colnames(S)))) > 0)
      stop(sprintf("%s in Q does not form part of the variable names of the data\n",
           Q[is.na(match(Q, colnames(S)))]))
    Q <- match(Q, colnames(S))
  }

  if (!R.code.only) {
    return(qpgraph:::.qpFastCItest(S, N, i, j, Q));
  }

  q       <- length(Q)
  Mmar    <- S[c(i, j, Q), c(i, j, Q)]
  S11     <- Mmar[1,1]
  S12     <- Mmar[1,-1]
  S21     <- Mmar[-1,1]
  S22     <- Mmar[-1,-1]
  S22inv  <- solve(S22)
  betahat <- S12 %*% S22inv[,1]
  sigma   <- sqrt((S11 - S12 %*% S22inv %*% S21) * (N - 1) / (N - q - 2))
  se      <- sigma * sqrt(S22inv[1,1] / (N - 1))
  t.value <- betahat / se
  p.value <- 2 * (1 - pt(abs(t.value), N - q - 2))

  return(list(t.value=t.value,p.value=p.value))
}



## function: qpHist
## purpose: plot the distribution of non-rejection rates
## parameters: nrrMatrix - matrix of non-rejection rates
##             K - concentration matrix from the generative graph
##             freq - logical; if TRUE, the histograms show frequencies (counts)
##                    of occurrence of the different non-rejection rate values;
##                    if FALSE, then probability densities are plotted
## return: none

qpHist <- function(nrrMatrix, K=NULL,
                   titlehist = "all estimated\nnon-rejection rates", freq=TRUE) {
  # all
  nrr <- nrrMatrix[upper.tri(nrrMatrix)]
  nrr_rg <- range(nrr)
  if(is.null(K)){
    hist(nrr, col="yellow", main=titlehist, xlab="non-rejection rate", freq=freq)
  } else {
    # only beta
    T <- nrrMatrix
    T[K == 0] <- -1
    xbeta <- T[upper.tri(T)]
    xbeta <- xbeta[xbeta != -1]
    # not beta
    T <- nrrMatrix
    T[K != 0] <- -1
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
    close.screen(all=TRUE)
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
##                           should be returned, either an incidence matrix,
##                           a matrix with the list of edges, a graphNEL structure
##                           or a graphAM structure
## return: incidence matrix of the qp-graph

qpGraph <- function(nrrMatrix, threshold=NULL, topPairs=NULL, pairup.i=NULL,
                    pairup.j=NULL, return.type=c("incidence.matrix", "edge.list",
                    "graphNEL", "graphAM")) {

  return.type <- match.arg(return.type)

  n.var <- nrow(nrrMatrix)

  if (is.null(colnames(nrrMatrix))) {
    vertex.labels <- as.character(1:n.var)
  } else {
    vertex.labels <- colnames(nrrMatrix)
  }

  if (is.null(threshold) && is.null(topPairs))
    stop("either threshold or topPairs must be set different to NULL\n")

  if (!is.null(threshold) && !is.null(topPairs))
    stop("only either threshold or topPairs can be set different to NULL\n")

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j must both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(nrrMatrix)))
      stop("when using pairup.i and pairup.j, nrrMatrix must have row and column names\n")

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

  # non-available NRRs imply no edges
  nrrMatrix[is.na(nrrMatrix)] <- 1.0

  if (!is.null(threshold)) {
    I <- nrrMatrix <= threshold
  } else { # topPairs
    nrrUppTriMatrix <- nrrMatrix[upper.tri(nrrMatrix)]
    rowUppTri <- row(nrrMatrix)[upper.tri(nrrMatrix)]
    colUppTri <- col(nrrMatrix)[upper.tri(nrrMatrix)]
    orderedMeasurementsIdx <- sort(nrrUppTriMatrix, index.return=TRUE,
                                   decreasing=FALSE)$ix
    ranking <- cbind(rowUppTri[orderedMeasurementsIdx],
                     colUppTri[orderedMeasurementsIdx])
    I <- matrix(FALSE, nrow=n.var, ncol=n.var)
    I[ranking[1:topPairs,]] <- TRUE
    I[cbind(ranking[1:topPairs,2], ranking[1:topPairs,1])] <- TRUE
  }

  rownames(I) <- colnames(I) <- vertex.labels
  diag(I) <- FALSE # whatever the threshold is the graph should have no loops

  if (return.type == "incidence.matrix") {
    return(I)
  } else if (return.type == "edge.list") {
    m <- cbind(row(I)[upper.tri(I) & I],col(I)[upper.tri(I) & I])
    colnames(m) <- c("i","j")
    return(m)
  } else if (return.type == "graphNEL") {
    require(graph)
    edL <- vector("list",length=n.var)
    names(edL) <- vertex.labels
    for (i in 1:n.var)
      edL[[i]] <- list(edges=(1:n.var)[I[i,]],weights=rep(1,sum(I[i,])))
    g <- new("graphNEL",nodes=vertex.labels,edgeL=edL,edgemode="undirected")
    return(g)
  } else if (return.type == "graphAM") {
    require(graph)
    g <- new("graphAM",adjMat=I+0,edgemode="undirected",values=list(weight=1))
    return(g)
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
##                           should be returned, either an incidence matrix,
##                           a matrix with the list of edges, a graphNEL structure
##                           or a graphAM structure
## return: incidence matrix of the qp-graph

qpAnyGraph <- function(measurementsMatrix, threshold=NULL, remove=c("below", "above"),
                       topPairs=NULL, decreasing=TRUE, pairup.i=NULL, pairup.j=NULL,
                       return.type=c("incidence.matrix", "edge.list", "graphNEL",
                                     "graphAM")) {

  remove <- match.arg(remove)
  return.type <- match.arg(return.type)

  n.var <- nrow(measurementsMatrix)

  if (is.null(colnames(measurementsMatrix))) {
    vertex.labels <- as.character(1:n.var)
  } else {
    vertex.labels <- colnames(measurementsMatrix)
  }

  if (is.null(threshold) && is.null(topPairs))
    stop("either threshold or topPairs must be set different to NULL\n")

  if (!is.null(threshold) && !is.null(topPairs))
    stop("only either threshold or topPairs can be set different to NULL\n")

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j must both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(measurementsMatrix)))
      stop("when using pairup.i and pairup.j, measurementsMatrix must have row and column names\n")

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

  # non-available measurements imply no edges
  measurementsMatrix[is.na(measurementsMatrix)] <- NA

  if (!is.null(threshold)) {
    if (remove == "below")
      I <- measurementsMatrix >= threshold
    else
      I <- measurementsMatrix <= threshold
  } else { # topPairs
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
    I <- matrix(FALSE, nrow=n.var, ncol=n.var)
    I[ranking[1:topPairs,]] <- TRUE
    I[cbind(ranking[1:topPairs,2], ranking[1:topPairs,1])] <- TRUE
  }

  I[is.na(I)] <- FALSE
  rownames(I) <- colnames(I) <- vertex.labels
  diag(I) <- FALSE # whatever the threshold is the graph should have no loops

  if (return.type == "incidence.matrix") {
    return(I)
  } else if (return.type == "edge.list") {
    m <- cbind(row(I)[upper.tri(I) & I],col(I)[upper.tri(I) & I])
    colnames(m) <- c("i","j")
    return(m)
  } else if (return.type == "graphNEL") {
    require(graph)
    edL <- vector("list",length=n.var)
    names(edL) <- vertex.labels
    for (i in 1:n.var)
      edL[[i]] <- list(edges=(1:n.var)[I[i,]],weights=rep(1,sum(I[i,])))
    g <- new("graphNEL",nodes=vertex.labels,edgeL=edL,edgemode="undirected")
    return(g)
  } else if (return.type == "graphAM") {
    require(graph)
    g <- new("graphAM",adjMat=I+0,edgemode="undirected",values=list(weight=1))
    return(g)
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
    colnames(matgdthr) <- c("density","threshold")
    n.var <- nrow(nrrMatrix)
    n.adj <- n.var*(n.var-1)/2

    for (i in 1:length(br)) {
      threshold <- br[i]
      nrrMatrix[is.na(nrrMatrix)] <- 1.0 # non-available NRRs imply no edges
      I <- nrrMatrix <= threshold
      diag(I) <- FALSE # if the threshold is 1.0 the resulting qp-graph
                       # will be the complete undirected graph but
                       # still it should have no loops
      n.edg <- sum(I) / 2
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

  return(list(data=matgdthr,sparseness=1-area$value))
}



## function: qpCliqueNumber
## purpose: calculate the size of the largest maximal clique in a given undirected graph
## parameters: g - either a graphNEL object or an incidence matrix of the graph
##             exact.calculation - flag that when set to TRUE the exact maximum clique
##                                 size is calculated and when set to FALSE a lower
##                                 bound is calculated instead
##             return.vertices - returns one set of vertices forming a maximal clique of the
##                               maximum size when this flag is set to TRUE
##             approx.iter - number of iterations performed to calculate
##                           the lower bound on the clique number of
##                           each graph (exact.calculation is FALSE)
##             verbose - show progress on the clique number calculation
## return: the size of the largest maximal clique in the given graph, also known as
##         its clique number

qpCliqueNumber <- function(g, exact.calculation=TRUE, return.vertices=FALSE,
                           approx.iter=100, verbose=TRUE) {

  if (class(g) == "graphNEL") {
    require(graph)
    if (edgemode(g) != "undirected")
      stop("g must be an undirected graph\n")

    I <- as(g, "matrix") == 1
  } else if (is.matrix(g)) {
    I <- g
    if (nrow(I) != ncol(I))
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!identical(I, t(I)))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g must be either a graphNEL object or a boolean incidence matrix\n")

  n.var <- nrow(I)
  n.possibleedges <- (n.var * (n.var-1)) / 2

  if (sum(I)/2 == 0) {
    return(1)
  }

  if (sum(I)/2 == n.possibleedges) {
    return(n.var)
  }

  maximum_clique <- 0

  if (exact.calculation == TRUE) {

    I <- I == 1 # make sure we get a boolean matrix

    maximum_clique <- qpgraph:::.qpCliqueNumberOstergard(I,return.vertices=return.vertices,verbose=verbose)
  } else {

    if (verbose) {
      cat("calculating lower bound on the maximum clique size\n")
    }

    clique.number <- 0
    clique.vertices <- c()

    I <- I + 0 # make sure we get a 0-1 matrix
    deg <- sort(rowsum(I,rep(1,n.var)),index.return=TRUE,decreasing=TRUE) # order by degree

    ppct <- -1
    for (i in 1:approx.iter) {

      pdeg <- deg$ix
      if (i %% n.var + 1 > 1) {
        sset <- sample(1:n.var,i %% n.var + 1,replace=FALSE) # we alter the order of the ranking
        ssetelem <- pdeg[sset]                               # by degree with increasing levels
        ssetelem <- sample(ssetelem)                         # of randomness cyclically
        pdeg[sset] <- ssetelem
      }
      clq <- c(pdeg[1])
      j <- 2
      for (j in 2:n.var) {
        v <- pdeg[j]
        clq2 <- c(clq,v)
        if (sum(I[clq2,clq2]) == length(clq2)*length(clq2)-length(clq2)) {
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

  return(maximum_clique)
}



## function: qpClique
## purpose: calculate and plot the maximum clique size as function of the
##          non-rejection rate
## parameters: nrrMatrix - matrix of non-rejection rates
##             N - sample size
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
##         sample size N, and the resulting maximum clique size

qpClique <- function(nrrMatrix, N=NA, threshold.lim=c(0,1), breaks=5, plot=TRUE,
                     exact.calculation=TRUE, approx.iter=100,
                     qpCliqueOutput=NULL, density.digits=0,
                     logscale.clqsize=FALSE,
                     titleclq="maximum clique size as function of threshold",
                     verbose=FALSE) {

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
    n.var <- nrow(nrrMatrix)
    n.adj <- n.var*(n.var-1)/2

    for (i in 1:length(br)) {
      if (verbose) {
        cat(paste("break: ",i,sep=""))
        cat("\n")
      }
      threshold <- br[i]
      nrrMatrix[is.na(nrrMatrix)] <- 1.0 # non-available NRRs imply no edges
      I <- nrrMatrix <= threshold
      diag(I) <- FALSE # if the threshold is 1.0 the resulting qp-graph
                       # will be the complete undirected graph but
                       # still it should have no loops
      n.edg <- sum(I) / 2
      dimnames(I) <- list(1:length(I[,1]), 1:length(I[1,]))
      maxsize <- qpCliqueNumber(I, exact.calculation, approx.iter=approx.iter,
                                  verbose=verbose)
      mpctedclqsze[i,] <- c((n.edg / n.adj) * 100, maxsize, threshold)
      if (!is.na(N)) {
        if (maxsize > maxclqszeunderN && maxsize < N) {
          maxclqszeunderN <- maxsize
          thrmaxclqszeunderN <- threshold
        }
      }
    }

    if (is.na(N))
      maxclqszeunderN <- thrmaxclqszunderN <- NA

  } else {
    br <- qpCliqueOutput$data[,3]
    mpctedclqsze <- qpCliqueOutput$data
    thrmaxclqszeunderN <- qpCliqueOutput$threshold
    maxclqszeunderN <- qpCliqueOutput$clqsizeunderN
    N <- qpCliqueOutput$N
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
         ylim=range(1,N,mpctedclqsze[,2],na.rm=TRUE),lwd=2,lty=linetype,
         xlab="threshold",ylab="maximum clique size",main=titleclq,log=logscale)
    if (!is.na(N))
      abline(h=N,col="red",lwd=2)
    pct <- round(mpctedclqsze[,1],digits=density.digits)
    text(br,mpctedclqsze[,2],lab=paste(pct,"%",sep=""),pos=1,cex=.7)
    legend(min(threshold.lim),max(N,mpctedclqsze[,2],na.rm=TRUE),label,lty=linetype,lwd=2,pch=1)
  }

  m <- mpctedclqsze[,c(3,2)]
  m[,2] <- m[,2] / n.var
  f <- approxfun(m)
  area <- integrate(f,min(m[,1]),max(m[,1]))

  return(list(data=mpctedclqsze,complexity=area$value,threshold=thrmaxclqszeunderN,
              clqsizeunderN=maxclqszeunderN,N=N,exact.calculation=exact.calculation))
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
## parameters: g - either a graphNEL object or an incidence matrix of the graph
##             clqspervtx - store the indices of the cliques where each vertex
##                          belongs to in the first p entries (|V|=p) of the
##                          returned list
##             verbose - show progress on the clique search
## return: a list of maximal cliques

qpGetCliques <- function(g, clqspervtx=FALSE, verbose=TRUE) {

  if (class(g) == "graphNEL") {
    require(graph)
    if (edgemode(g) != "undirected")
      stop("g must be an undirected graph\n")

    I <- as(g, "matrix") == 1
  } else if (is.matrix(g)) {
    I <- g
    if (nrow(I) != ncol(I))
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!identical(I, t(I)))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g must be either a graphNEL object or a boolean incidence matrix\n")

  return(qpgraph:::.qpFastCliquerGetCliques(I,clqspervtx=clqspervtx,verbose=verbose))
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

  if (!R.code.only) {
    return(qpgraph:::.qpFastIPF(vv,clqlst,tol,verbose))
  }

  if (verbose) {
    n.var <- nrow(vv)
    if (clqlst[[1]][1] <= n.var) {
      n.var <- 0
    }
    cat(paste(paste("qpIPF: ",length(clqlst)-n.var),"cliques\n"))
  }

  V <- diag(length(vv[, 1]))
  Vold <-  - V
  while(max(abs(V - Vold)) > tol) {
    Vold <- V
    V <- qpgraph:::.qpIPFpass(vv, V, clqlst)
    if (verbose)
      cat("qpIPF: precision =", max(abs(V - Vold)), "\n")
  }

  return(V)
}



## function: qpPAC
## purpose: for a given undirected graph in an incidence matrix estimate the
##          partial correlation coefficient (PAC) and its corresponding p-value
##          for each edge in the graph
## parameters: data - data set from where to estimate the PACs
##             g - either a graphNEL object or an incidence matrix of the graph
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             return.K - flag that when set to TRUE the function also returns
##                        the concentration matrix K; if FALSE (default) does not
##                        return K
##             verbose - flag that when set to TRUE the IPF algorithm
##                       shows the convergence progression
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two matrices, one with the estimates of the PACs and
##         the other with their p-values

setGeneric("qpPAC", function(data, ...) standardGeneric("qpPAC"))

# data comes as an ExpressionSet object
setMethod("qpPAC", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpPAC(S, N, ...)
          })

# data comes as a data frame
setMethod("qpPAC", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpPAC(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpPAC", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            S <- cov(data)
            N <- length(data[,1])
            if (!is.null(colnames(data))) 
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpPAC(S, N, ...)
          })

.qpPAC <- function(S, N, g, return.K=FALSE, verbose=TRUE, R.code.only=FALSE) {

  if (class(g) == "graphNEL") {
    require(graph)
    if (edgemode(g) != "undirected")
      stop("g must be an undirected graph\n")

    I <- as(g, "matrix") == 1
  } else if (is.matrix(g)) {
    I <- g
    if (nrow(I) != ncol(I))
      stop("g is not an squared matrix nor a graphNEL object\n")

    if (!identical(I, t(I)))
      stop("g is not a symmetric matrix nor a graphNEL object\n")
  } else
    stop("g must be either a graphNEL object or a boolean incidence matrix\n")

  var.names <- rownames(S)
  n.var <- nrow(S)
  dimnames(I) <- dimnames(S)

  # get the cliques

  clqlst <- qpGetCliques(I,verbose=verbose)

  # get a maximum likelihood estimate of the sample covariance matrix
  # using the clique list

  Shat <- qpIPF(S,clqlst,verbose=verbose,R.code.only=R.code.only)

  # estimate partial correlation coefficients and their standard errors

  K <- solve(Shat)
  SE <- qpgraph:::.qpEdgePACSE(Shat, I, R.code.only=R.code.only)

  # return matrices of partial correlations, standard errors
  # and p-values for every edge

  C <- N * (K^2 / SE)
  offdiag_minsgn_mask <- -1 * (!diag(n.var)) + diag(n.var)
  rho_coef <- offdiag_minsgn_mask * qpDscale(K)
  p.values <- 1 - pchisq(C, df=1)
  dimnames(rho_coef) <- dimnames(p.values) <- list(var.names, var.names)

  list2return <- list(R=rho_coef, P=p.values)
  if (return.K)
    list2return <- list(R=rho_coef, P=p.values, K=K)

  return(list2return)
}



## function: qpPCC
## purpose: estimate pairwise Pearson correlation coefficients (PCCs) between all
##         pairs of variables
## parameters: data - data set from where to estimate the PCCs
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
## return: a list with two matrices, one with the estimated PCCs and the
##         other with their p-values

setGeneric("qpPCC", function(data, ...) standardGeneric("qpPCC"))

# data comes as an ExpressionSet object
setMethod("qpPCC", signature(data="ExpressionSet"),
          function(data, ...) {
            exp <- t(exprs(data))
            S <- cov(exp)
            N <- length(sampleNames(data))
            rownames(S) <- colnames(S) <- featureNames(data)
            qpgraph:::.qpPCC(S, N, ...)
          })

# data comes as a data frame
setMethod("qpPCC", signature(data="data.frame"),
          function(data, long.dim.are.variables=TRUE, ...) {
            m <- as.matrix(data)
            rownames(m) <- rownames(data)
            if (!is.double(m))
              stop("data should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(m),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              m <- t(m)
            S <- cov(m)
            N <- length(m[,1])
            if (!is.null(colnames(m)))
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpPCC(S, N, ...)
          })

          
# data comes as a matrix
setMethod("qpPCC", signature(data="matrix"),
          function(data, long.dim.are.variables=TRUE, ...) {
            if (long.dim.are.variables &&
              sort(dim(data),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              data <- t(data)

            S <- cov(data)
            N <- length(data[,1])
            if (!is.null(colnames(data))) 
              rownames(S) <- colnames(S) <- 1:nrow(S)
            else
              rownames(S) <- colnames(S) <- colnames(data)
            qpgraph:::.qpPCC(S, N, ...)
          })

.qpPCC <- function(S, N) {

  var.names <- rownames(S)

  # estimate PCCs by scaling the covariance matrix
  R <- qpDscale(S)

  # calculate t-statistics
  T <- (N - 2) / (1 - R*R)
  diag(T) <- (N - 2) * 100000 # just to get 0 p-values on the diagonal
  T <- R * sqrt(T)

  # calculate two-sided p-values
  p <- pt(T, df=N - 2)
  P <- 2 * pmin(p, 1 - p)

  list(R=R,P=P)
}



## function: qpRndGraph
## purpose: builds a random undirected graph with n.vtx vertices
##          and for every vertex its boundary <= n.bd
## parameters: n.vtx - number of vertices
##             n.bd - maximum boundary for every vertex
## return: the incidence matrix of the resulting graph

qpRndGraph <- function(n.vtx, n.bd){
  #
  # function to remove the zeros at the end of a vector.
  #
  ext.vec <- function(v){v=v[0:sum(v!=0)]}
  G <- matrix(0, n.vtx, n.bd)
  vtx.list <- 1:n.vtx
  full.bd  <- c() 
  for(i in 1:n.vtx){
    #
    # computes n.sam the number of vertices to add to 
    # the boundary of i and the the population vtx.pop 
    # of candidate vertices to be in the boundary of i
    #
    bd.i  <- ext.vec(G[i,])
    n.sam <- n.bd-length(bd.i)
    vtx.pop <- vtx.list[-c(i, bd.i, full.bd)]
    #
    # update the boundary. Note that sample() works 
    # in a different way if l.vp==1
    #
    l.vp=length(vtx.pop)
    if(l.vp>1 & l.vp>=n.sam){
        new.bd  <- sample(x=vtx.pop, size=n.sam, rep=F) 
        G[i,] <- c(bd.i, new.bd)
    } else {
        if(n.sam>0) {new.bd<- vtx.pop}else{new.bd<-c()}
        G[i,] <- c(bd.i, new.bd, rep(0, length=(n.bd-length(c(bd.i, new.bd)))))
    }
    #
    # updates the list of vertices with full boundary
    #
    full.bd <- c(full.bd, i)
    #
    # add the vertex i to the boundary of vertex added to the boundary of i
    #
    for (j in new.bd){
      e <- ext.vec(G[j,])
      if (length(e)==(n.bd-1)){
        full.bd <- c(full.bd, j)
        G[j, ] <- c(e, i)
      }else{
        G[j, ] <- c(e, i, rep(0, length=(n.bd-length(e)-1)) )
      }

    }
  }
  #
  # denotes by n.vtx+1 all the non-assigned boundary locations.  
  #
  G[G==0] <- n.vtx+1

  # get the corresponding incidence matrix
  I <- matrix(FALSE, nrow=n.vtx, ncol=n.vtx)

  for (i in 1:n.vtx){
      for(j in 1:n.bd){
          if(G[i, j]!=(n.vtx+1)){
              I[i, G[i, j]] <- I[G[i, j], i]<- TRUE
          }             
      }
  }

  return(I)
}



## function: qpSampleMvnorm
## purpose: sample N independent observations from a multivariate normal
##          distribution with a given mean vector and concentration matrix
## paramters: K - concentration matrix
##            N - number of observations to sample
##            mean - mean vector
## return: a matrix where rows correspond to observations and columns
##         to random variables

qpSampleMvnorm <- function(K, N, mean = rep(0, nrow(K))) {
  require(mvtnorm)

  n.var <- dim(K)[1]
  sigma <- qpDscale(solve(K))
  X <- rmvnorm(N, mean, sigma)

  return(X)
}



## function: qpI2K
## purpose: builds a random concentration matrix from an incidence matrix
## parameters: I - incidence matrix
##             verbose - output progress
##             R.code.only - flag set to FALSE when using the C implementation
## return: a random concentration matrix with zeroes at the empty adjacencies of
##         the undirected graph defined by the input incidence matrix I

qpI2K <- function(I, verbose=FALSE, R.code.only=FALSE) {
  n.var <- nrow(I)

  right <- FALSE
  while (!right) {
    right <- TRUE

    ## generate a random correlation matrix
    ## the four lines below were adapted from the code
    ## for the rcorr function from the ggm package by
    ## G. Marchetti implementing the method from Marsaglia & Oltkin
    K <- matrix(rnorm(n.var*n.var),nrow=n.var,ncol=n.var)
    d <- apply(K,1,function(x) sqrt(sum(x*x)))
    K <- sweep(K,1,d,"/")
    K <- K %*% t(K)

    S <- qpDscale(solve(K))

    clqlst <- qpGetCliques(I,verbose=verbose)
    Shat <- qpIPF(S,clqlst,verbose=verbose,R.code.only=R.code.only)
    Khat <- qpDscale(solve(Shat))
    Khat[abs(Khat) < 0.001] <- 0

    if(sum(Khat[!I]!=0)!= n.var || sum(diag(solve(Khat)) < 0) > 0) {
      warning("something wrong in the zero structure of K, trying again")
      right <- FALSE
    }
  }

  return(Khat)
}



## function: qpK2R
## purpose: obtain the partial correlation coefficients from a given
##          concentration matrix
## parameters: K - concentration matrix
## return: a matrix with the partial correlation coefficients

qpK2R <- function(K) {
  n.var <- nrow(K)
  offdiag_minsgn_mask <- -1 * (!diag(n.var)) + diag(n.var)
  R <- offdiag_minsgn_mask * qpDscale(K)
  return(R)
}



## function: qpPrecisionRecall
## purpose: calculate the precision-recall curve for a given measure of
##          association between all pairs of variables in a matrix
## parameters: measurementsMatrix - matrix containing the measure of association
##                                  between all pairs of variables
##             refI - incidence matrix of reference from which to calculate
##                    the precision-recall curve
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

qpPrecisionRecall <- function(measurementsMatrix, refI, decreasing=TRUE,
                              pairup.i=NULL, pairup.j=NULL,
                              recallSteps=c(seq(0,0.1,0.005),seq(0.2,1.0,0.1))) {

  if (!is.matrix(measurementsMatrix))
    stop("measurementsMatrix must be a matrix\n")

  if (!is.matrix(refI))
    stop("refI must be a matrix\n")

  if (nrow(measurementsMatrix) != ncol(measurementsMatrix))
    stop("measurementsMatrix must be a squared matrix\n")

  if (nrow(refI) != ncol(refI))
    stop("refI must be a squared matrix\n")

  if (nrow(measurementsMatrix) != nrow(refI) ||
      ncol(measurementsMatrix) != ncol(refI))
    stop("measurementsMatrix and refI must have the same dimensions\n")

  n.var <- nrow(measurementsMatrix)

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j must both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(measurementsMatrix)))
      stop("when using pairup.i and pairup.j, measurementsMatrix must have row and column names\n")

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
  total_positives <- sum(refI[upper.tri(refI) & !is.na(measurementsMatrix)])

  status <- refI[as.matrix(edgeRnk[,c(1,2)])]
  status_tp <- rep(0, length(status))
  status_tp[status] <- 1:total_positives
  preRec <- matrix(0, nrow=length(recallSteps), ncol=4)
  colnames(preRec) <- c("Recall","Precision","TP","ScoreThreshold")

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

    preRec[i, ] <- c(actualRecall, precision, tp, scoreThreshold)
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



## function: qpDscale
## purpose: scale a matrix by its diagonal
## parameters: V - the matrix
## return: the scaled matrix

qpDscale <- function(V){
  d <- 1 / sqrt(diag(V))
                                                                                                  
  return(V * outer(d, d))
}



## function: qpFunctionalCoherence
## purpose: estimate functional coherence of a transcripcional regulatory network
##          represented by means of an undirected graph encoded by an incidence
##          matrix and of a set of transcription factor genes. In these
##          calculations Gene Ontology (GO) annotations are employed through a
##          given annotation .db package for the Entrez Gene IDs associated to
##          the rows and columns of the incidence matrix.
## parameters: I - incidence matrix of the undirected graph representing the
##                 transcriptional regulatory network
##             TFgenes - vector of transcription factor gene names (matching the
##                       genes at the rows and column names of I)
##             chip - name of the .db package containing the GO annotations
##             minRMsize - minimum size of the target gene set in each regulatory
##                         module where functional enrichment will be calculated
##                         and thus where functional coherence will be estimated
##             verbose - logical; if TRUE the function will show progress on the
##                       calculations; if FALSE will remain quiet (default)
## return: a list with three slots, a first one containing the transcriptional
##         regulatory network as a list of regulatory modules and their targets,
##         a second one containing this same network but including only those
##         modules with GO BP annotations and a third one consisting of a vector
##         of functional coherence values

qpFunctionalCoherence <- function(I, TFgenes, chip, minRMsize=5, verbose=FALSE) {
  require(GOstats)

  if (is.null(colnames(I)) || is.null(rownames(I)))
    stop("incidence matrix I must have row and column names corresponding to the gene IDs")

  if (length(TFgenes) < 1)
    stop("TFgenes must contain at least one transcription factor gene\n")

  allGenes <- rownames(I)

  if (!is.character(TFgenes))
    stop("gene identifiers in TFgenes must belong to the class character\n")

  if (sum(is.na(match(TFgenes, allGenes))) > 0)
    stop("TFgenes is not a subset from the genes defining the incidence matrix I\n")

  p <- dim(I)[1]
  geneBPuniverse <- qpgraph:::.qpFilterByGO(allGenes, chip, "BP")

  TFgenes_i <- match(TFgenes, allGenes)
  txRegNet <- lapply(TFgenes_i, function(x) allGenes[I[as.integer(x), ]])
  names(txRegNet) <- TFgenes
  regModuleSize <- unlist(lapply(txRegNet, length))
  names(regModuleSize) <- TFgenes

  if (verbose)
    cat(sprintf("qpFunctionalCoherence: calculating GO enrichment in %d target gene sets\n",
        length(TFgenes[regModuleSize >= minRMsize])))

  txRegNetGO <- list()
  for (TFgene in TFgenes[regModuleSize >= minRMsize]) {
    if (verbose)
      cat(".")
    TFgeneTGs <- txRegNet[[TFgene]]
    TFgeneTGsWithGO <- qpgraph:::.qpFilterByGO(TFgeneTGs, chip, "BP")

    if (length(TFgeneTGsWithGO) >= minRMsize) {
      goHypGparams <- new("GOHyperGParams",
                          geneIds=TFgeneTGsWithGO,
                          universeGeneIds=geneBPuniverse,
                          annotation=chip, ontology="BP",
                          pvalueCutoff=0.05, conditional=TRUE,
                          testDirection="over")
      goHypGcond <- hyperGTest(goHypGparams)
      txRegNetGO[[TFgene]] <- list(initialRMsize=regModuleSize[TFgene],
                                   TGgenesWithGO=TFgeneTGsWithGO,
                                   goBPcondResult=goHypGcond,
                                   goBPcondResultSigCat=sigCategories(goHypGcond))
    }
  }

  if (verbose)
    cat(sprintf("\nqpFunctionalCoherence: calculating functional coherence in %d RMs\n",
        length(names(txRegNetGO))))

  TFgenesWithGO <- qpgraph:::.qpFilterByGO(TFgenes, chip, "BP")
  TFgenesWithGO <- AnnotationDbi::mget(TFgenesWithGO,get(gsub(".db","GO",chip)))
  TFgenesWithGO <- lapply(TFgenesWithGO,
                          function(x) if (is.list(x)) {
                                        z <- sapply(x, function(x) x$Ontology);
                                        z[unique(names(z))]
                                      })
  TFgenesWithGOBP <- lapply(TFgenesWithGO,
                            function(x) if (sum(x=="BP",na.rm=TRUE) > 0) {
                                          names(x[x=="BP" & !is.na(x)])
                                        } else { NULL })

  goTerms <- unlist(AnnotationDbi::eapply(GOTERM, function(x) x@Term))
  goTermOntologies <- unlist(AnnotationDbi::eapply(GOTERM, function(x) x@Ontology))
  goTermBPOntology <- names(goTermOntologies[goTermOntologies == "BP"])

  # remove from the GO annotation of the transcription factor those GO terms
  # that have the word 'transcription', i.e., we try to remove terms associated
  # to transcriptional regulation from the transcription factor GO annotation
  TFgoTerms <- names(goTerms[grep("transcription", goTerms)])
  TFgoTerms <- TFgoTerms[!is.na(match(TFgoTerms,goTermBPOntology))]

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
        if (length(txRegNetGO[[TFgene]]$goBPcondResultSigCat) == 0)
          sUI <- 0
        else { # otherwise the functional coherence is estimated as the similarity
               # between the GO graphs associated to the transcription factor GO
               # annotations and the GO over-represented terms in the target gene set
          gTF <- GOGraph(TFgoAnnot, GOBPPARENTS)
          gTF <- removeNode("all", gTF)
          gTG <- GOGraph(txRegNetGO[[TFgene]]$goBPcondResultSigCat, GOBPPARENTS)
          gTG <- removeNode("all", gTG)
          sUI <- simUI(gTF, gTG)
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
##             I - incidence matrix of the graph and thus it is assumed that the diagonal
##                 is set to either 0s or FALSE truth values since there should be no loops
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two members: K - the concentration matrix; SE the matrix
##         with the standard errors of the edges

.qpEdgePACSE <- function(S, I, R.code.only=FALSE) {

  if (!R.code.only) {
    return(qpgraph:::.qpFastPACSE(S, I));
  }

  n.var <- nrow(I)

  I <- I + diag(n.var) # in the code below we need 1s in the main diagonal and
                       # then at the same time we make sure we get a 0-1 matrix
                       # as a truth value + 0 or 1 equals a number

  I[col(I) > row(I)] <- NA

  # selection row and column indices corr. to the non-zero elem.
  I[I == 0] <- NA
  r.nz <- c(row(I))[!is.na(I)]
  c.nz <- c(col(I))[!is.na(I)]

  # computation of the Fisher information matrix
  Iss <- S[c.nz,c.nz] * S[r.nz,r.nz] +
         S[c.nz,r.nz] * S[r.nz,c.nz]
  Iss <- solve(Iss)
  IssI <- matrix(rep(0, length(r.nz) * length(c.nz)), nrow=length(r.nz))
  diag(IssI) <- 1
  IssI[cbind((1:length(r.nz))[r.nz==c.nz], (1:length(r.nz))[r.nz==c.nz])] <- 2

  FISHER <- IssI %*% Iss %*% IssI

  # standard errors are in the diagonal of the Fisher information matrix
  F <- diag(FISHER)
  SE <- matrix(NA, nrow(I), nrow(I))
  SE[cbind(r.nz,c.nz)] <- SE[cbind(c.nz,r.nz)] <- F
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

  haveGo <- sapply(AnnotationDbi::mget(entrezGeneIds,
                                       getAnnMap(map="GO", chip=chip, type="db")),
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



########################################################################
# PRIVATE FUNCTIONS THAT ARE ENTRY POINTS OF THE C CODE OF THE PACKAGE #
########################################################################

.qpFastNrr <- function(S, N, q, nTests, alpha, pairup.i.noint, pairup.j.noint,
                       pairup.ij.int, verbose) {
  return(.Call("qp_fast_nrr",S,as.integer(N),as.integer(q),as.integer(nTests),
                             as.double(alpha),as.integer(pairup.i.noint),
                             as.integer(pairup.j.noint),as.integer(pairup.ij.int),
                             as.integer(verbose)))
}

.qpFastEdgeNrr <- function(S, N, i, j, q, nTests, alpha) {
  return(.Call("qp_fast_edge_nrr",S,as.integer(N),as.integer(i),as.integer(j),
                                  as.integer(q),as.integer(nTests),
                                  as.double(alpha)))
}

.qpFastCItest <- function(S, N, i, j, C=c()) {
  return(.Call("qp_fast_ci_test",S,as.integer(N),as.integer(i),as.integer(j),C))
}

.qpFastCliquerGetCliques <- function(I,clqspervtx,verbose) {
  return(.Call("qp_fast_cliquer_get_cliques",I,clqspervtx,verbose))
}

.qpFastPACSE <- function(Shat, I) {
  return(.Call("qp_fast_pac_se",Shat ,I))
}

.qpFastIPF <- function(vv, clqlst, tol = 0.001, verbose = FALSE) {
  return(.Call("qp_fast_ipf",vv,clqlst,tol,verbose))
}

.qpCliqueNumberOstergard <- function(I,return.vertices,verbose) {
 return(.Call("qp_clique_number_os",I,return.vertices,verbose))
}
