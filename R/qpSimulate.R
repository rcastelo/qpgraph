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



## function: qpRndGraph
## purpose: samples a d-regular graph uniformly at random
## parameters: p - number of vertices
##             d - maximum boundary for every vertex
##             labels - vertex labels
##             exclude - vector of vertices inducing edges that should be
##                       excluded from the sampled d-regular graph
##             verbose - show progress on the simulation
##             return.type - type of object to return
## return: the adjacency matrix of the resulting graph

qpRndGraph <- function(p=6, d=2, labels=1:p, exclude=NULL, verbose=FALSE,
                       return.type=c("adjacency.matrix", "edge.list", "graphBAM", "graphNEL"), R.code.only=FALSE) {
  return.type <- match.arg(return.type)

  if ((p*d) %% 2 != 0)
    stop("The number of vertices p times the degree d of each vertex, i.e., the product p x d, should be even in order to sample a d-regular graph on p vertices uniformly at random\n")

  if (d > sqrt(p))
    warning("Steger and Wormald (1999, see help page of this function) believe that when d >> sqrt(p) the resulting d-regular graph on p vertices may no longer be sampled from a approximately uniform probability distribution.\n")

  if (!is.null(exclude)) {
    if (!is.integer(exclude) || any(exclude < 1))
      stop("The argument exclude must be a vector of positive integers\n")
    if (any(is.na(match(exclude, 1:p))))
      stop("The argument exclude contains vertices outside range 1:p\n")
    if (any(duplicated(exclude)))
      stop("The argument exclude contains duplicated vertices\n")
  }

  G <- matrix(FALSE, nrow=p, ncol=p)

  if (!R.code.only) {
    G <- .qpFastRndGraph(p, d, exclude, verbose)
  } else {

    if (verbose)
      pb <- txtProgressBar(style=3, min=0, max=p)

    Gex <- matrix(FALSE, nrow=p, ncol=p)
    if (!is.null(exclude)) {
      Gex[exclude, exclude] <- TRUE
    }
  
    while (any(rowSums(G) != d)) {
      G <- matrix(FALSE, nrow=p, ncol=p, dimnames=list(1:p, 1:p))
      S <- TRUE
    
      while (!all(is.na(S))) {
        dG <- rowSums(G) ## calculate degree
        dG[dG > d-1] <- NA ## select pairs where both vertices have
                           ## degree at most d-1
        S <- (d - dG) %o% (d - dG)
        S[G] <- NA ## exclude adjacent pairs of vertices

        if (!all(is.na(S))) { ## if there are missing edges to add
          S <- S / sum(S[upper.tri(S) & !Gex], na.rm=TRUE)
          ridx <- row(S)[upper.tri(S) & !is.na(S) & !Gex]
          cidx <- col(S)[upper.tri(S) & !is.na(S) & !Gex]
          S <- S[upper.tri(S) & !is.na(S) & !Gex]
          cdf <- sort(S, decreasing=TRUE, index.return=TRUE) ## build CDF
          r <- runif(1, min=0, max=1)
          i <- cdf$ix[sum(r > cumsum(cdf$x))+1] ## sample one edge from the CDF
          G[ridx[i], cidx[i]] <- G[cidx[i], ridx[i]] <- TRUE ## add it
        }
  
        if (verbose)
          setTxtProgressBar(pb, sum(rowSums(G) == d))
      }
    }

    if (verbose)
      close(pb)
  }

  if (return.type == "edge.list") {
    m <- cbind(labels[row(G)[upper.tri(G) & G]], labels[col(G)[upper.tri(G) & G]])
    colnames(m) <- c("i", "j")
    return (m)
  } else if (return.type == "graphBAM") {
    df <- data.frame(from=labels[row(G)[upper.tri(G) & G]],
                     to=labels[col(G)[upper.tri(G) & G]],
                     weight = rep(1, sum(G)/2))
    Gbam <- graphBAM(df, edgemode = "undirected", nodes=labels)
    return (Gbam)
  } else if (return.type == "graphNEL") {
    df <- data.frame(from=labels[row(G)[upper.tri(G) & G]],
                     to=labels[col(G)[upper.tri(G) & G]],
                     weight = rep(1, sum(G)/2))
    Gnel <- as(graphBAM(df, edgemode = "undirected", nodes=labels), "graphNEL")
    return (Gnel)
  }


  ## return.type == "adjacency.matrix" using a memory-efficient lspMatrix-classs object
  G <- as(G, "lspMatrix")
  dimnames(G) <- list(1:p, 1:p)

  return(G)
}


.qpRndRegularGraph <- function(p=6, d=2, labels=1:p, exclude=NULL, verbose=FALSE,
                               return.type=c("adjacency.matrix", "edge.list", "graphBAM", "graphNEL"), R.code.only=FALSE) {
  return.type <- match.arg(return.type)

  if ((p*d) %% 2 != 0)
    stop("The number of vertices p times the degree d of each vertex, i.e., the product p x d, should be even in order to sample a d-regular graph on p vertices uniformly at random\n")

  if (d > sqrt(p))
    warning("Steger and Wormald (1999, see help page of this function) believe that when d >> sqrt(p) the resulting d-regular graph on p vertices may no longer be sampled from a approximately uniform probability distribution.\n")

  if (!is.null(exclude)) {
    if (!is.integer(exclude) || any(exclude < 1))
      stop("The argument exclude must be a vector of positive integers\n")
    if (any(is.na(match(exclude, 1:p))))
      stop("The argument exclude contains vertices outside range 1:p\n")
    if (any(duplicated(exclude)))
      stop("The argument exclude contains duplicated vertices\n")
  }

  G <- matrix(FALSE, nrow=p, ncol=p)

  if (!R.code.only) {
    G <- .qpFastRndGraph(p, d, exclude, verbose)
  } else {

    if (verbose)
      pb <- txtProgressBar(style=3, min=0, max=p)

    Gex <- matrix(FALSE, nrow=p, ncol=p)
    if (!is.null(exclude)) {
      Gex[exclude, exclude] <- TRUE
    }
  
    while (any(rowSums(G) != d)) {
      G <- matrix(FALSE, nrow=p, ncol=p, dimnames=list(1:p, 1:p))
      S <- TRUE
    
      while (!all(is.na(S))) {
        dG <- rowSums(G) ## calculate degree
        dG[dG > d-1] <- NA ## select pairs where both vertices have
                           ## degree at most d-1
        S <- (d - dG) %o% (d - dG)
        S[G] <- NA ## exclude adjacent pairs of vertices

        if (!all(is.na(S))) { ## if there are missing edges to add
          S <- S / sum(S[upper.tri(S) & !Gex], na.rm=TRUE)
          ridx <- row(S)[upper.tri(S) & !is.na(S) & !Gex]
          cidx <- col(S)[upper.tri(S) & !is.na(S) & !Gex]
          S <- S[upper.tri(S) & !is.na(S) & !Gex]
          cdf <- sort(S, decreasing=TRUE, index.return=TRUE) ## build CDF
          r <- runif(1, min=0, max=1)
          i <- cdf$ix[sum(r > cumsum(cdf$x))+1] ## sample one edge from the CDF
          G[ridx[i], cidx[i]] <- G[cidx[i], ridx[i]] <- TRUE ## add it
        }
  
        if (verbose)
          setTxtProgressBar(pb, sum(rowSums(G) == d))
      }
    }

    if (verbose)
      close(pb)
  }

  if (return.type == "edge.list") {
    m <- cbind(labels[row(G)[upper.tri(G) & G]], labels[col(G)[upper.tri(G) & G]])
    colnames(m) <- c("i", "j")
    return (m)
  } else if (return.type == "graphBAM") {
    df <- data.frame(from=labels[row(G)[upper.tri(G) & G]],
                     to=labels[col(G)[upper.tri(G) & G]],
                     weight = rep(1, sum(G)/2))
    Gbam <- graphBAM(df, edgemode = "undirected", nodes=labels)
    return (Gbam)
  } else if (return.type == "graphNEL") {
    df <- data.frame(from=labels[row(G)[upper.tri(G) & G]],
                     to=labels[col(G)[upper.tri(G) & G]],
                     weight = rep(1, sum(G)/2))
    Gnel <- as(graphBAM(df, edgemode = "undirected", nodes=labels), "graphNEL")
    return (Gnel)
  }


  ## return.type == "adjacency.matrix" using a memory-efficient lspMatrix-classs object
  G <- as(G, "lspMatrix")
  dimnames(G) <- list(1:p, 1:p)

  return(G)
}



.qpFastRndGraph <- function(p, d, exclude, verbose) {
  ## return(new("lspMatrix", Dim=c(as.integer(p), as.integer(p)),
  ##            Dimnames=list(1:p, 1:p),
  ##            x = .Call("qp_fast_rnd_graph", as.integer(p), as.integer(d),
  ##                      as.integer(exclude), as.integer(verbose))))
  A <- matrix(FALSE, nrow=p, ncol=p, dimnames=list(1:p, 1:p))
  A[upper.tri(A, diag=TRUE)] <- .Call("qp_fast_rnd_graph", as.integer(p), as.integer(d),
                                      as.integer(exclude), as.integer(verbose))
  A <- A | t(A)
  return(A)
}



## function: qpUnifRndAssociation
## purpose: builds a matrix of uniformly random correlation values between -1 and +1
## parameters: n.var - number of variables
##             var.names - names of the variables to put as dimension names
## return: a matrix of uniformly random correlation values between -1 and +1 for
##         every pair of variables

qpUnifRndAssociation <- function (n.var, var.names=as.character(1:n.var)) {
  n.var <- as.integer(n.var)
  x=runif((n.var*(n.var-1))/2+n.var, min=-1, max=+1)
  x[cumsum(1:n.var)] <- 1
  rndcor <- new("dspMatrix", Dim=c(n.var, n.var),
                Dimnames=list(var.names, var.names), x=x)

  return(rndcor)
}



## function: qpRndWishart
## purpose: Random generation for the (n.var x n.var) Wishart distribution with
##          matrix parameter A=diag(delta)%*%P%*%diag(delta) and degrees of
##          freedom df
## parameters: delta - a numeric vector of n.var positive values. If a scalar
##                     is provided then this is extended to form a vector.
##             P - a (n.var x n.var) positive definite matrix with unit
##                 diagonal. If a scalar is provided then this number is used
##                 as constant off-diagonal entry for P
##             df - degrees of freedom
##             n.var - dimension of the Wishart matrix. It is required only when
##                     both delta and P are scalar
## return: a list of two (n.var x n.var) matrices rW and meanW where rW is a
##         random value from the Wishart and meanW is the expected value of the
##         distribution

qpRndWishart <- function(delta=1, P=0, df=NULL, n.var=NULL) {

  if (length(delta) == 1 && length(P) == 1 && is.null(n.var))
    stop("The value of n.var is not specified and both delta and P are scalar")

  if (is.null(n.var)) n.var=max(length(delta), dim(P)[1])

  if (length(P) == 1) {
    P <- matrix(P, n.var, n.var)
    diag(P) <- 1
  } 

  if (max(abs(P) > 1) || min(eigen(P)$values)<=0 || !identical(P, t(P))) {
    stop("P should be either a (symmetric) positive definite matrix\n or a scalar larger than -(n.var-1)^(-1) and smaller than 1")
  } 

  if (length(delta) == 1) delta=rep(delta, length=n.var)
  if (min(delta) <= 0) stop("All entries of delta should be positive")
  if (is.null(df)) df <- n.var
  if (df <= (n.var-1)) stop("The value of df should be larger than (n.var-1)")

  Delta <- diag(delta)
  V <- Delta %*% P %*% Delta
  CV <- chol(V)
  CWS <- matrix(0, n.var, n.var)
  CWS[row(CWS) < col(CWS)] <- rnorm(n.var * (n.var - 1) / 2)
  diag(CWS) <- sqrt(rchisq(n.var, (df + 1)-(1:n.var)))
  CW <- CWS %*% CV
  W <- t(CW) %*% CW

  return(list(rW=as(W, "dspMatrix"), meanW=df * V))
}



## function: qpG2Sigma
## purpose: builds a random covariance matrix from an undirected graph
## parameters: g - undirected graph (either adjacency matrix, graphNEL object, graphAM object or graphBAM object)
##             rho - real number between 1/(n.var-1) and 1
##             matrix.completion - algorithm to perform matrix completion operations
##             verbose - output progress
##             R.code.only - flag set to FALSE when using the C implementation
## return: a random covariance matrix whose inverse contains zeroes at the
##         missing edges in G

qpG2Sigma <- function (g, rho=0, matrix.completion=c("HTF", "IPF"), tol=0.001,
                       verbose=FALSE, R.code.only = FALSE) {
  matrix.completion <- match.arg(matrix.completion)
  n.var <- NULL
  var.names <- NULL
  if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
    n.var <- nrow(g)
    var.names <- rownames(g)
    if (is.null(var.names))
      var.names <- 1:n.var
  } else if (class(g) == "graphNEL" || class(g) == "graphAM" || class(g) == "graphBAM") {
    n.var <- length(graph::nodes(g))
    var.names <- nodes(g)
  }

  if (is.null(n.var))
    stop("'g' is neither an adjacency matrix, graphNEL, graphAM or graphBAM type of object.\n")

  W <- qpRndWishart(delta=sqrt(1 / n.var), P=rho, n.var=n.var)$rW

  Sigma <- NULL
  if (matrix.completion == "IPF") {
    clqlst <- qpGetCliques(g, verbose=verbose)
    Sigma <- qpIPF(W, clqlst, tol=tol, verbose=verbose, R.code.only=R.code.only)
  } else
    Sigma <- qpHTF(W, g, tol=tol, verbose=verbose, R.code.only=R.code.only)

  Sigma <- as(Sigma, "dspMatrix")
  rownames(Sigma) <- colnames(Sigma) <- var.names

  return(Sigma)
}

