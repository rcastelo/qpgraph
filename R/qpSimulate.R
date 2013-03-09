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
    G <- qpgraph:::.qpFastRndGraph(p, d, exclude, verbose)
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


qpRndRegularGraph <- function(p=6, d=2, labels=1:p, exclude=NULL, verbose=FALSE,
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
    G <- qpgraph:::.qpFastRndGraph(p, d, exclude, verbose)
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

qpUnifRndAssociation <- function (n.var, var.names=1:n.var) {
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

oldrUGgmm <- function(n=1L, p=5L, d=2L, labels=sprintf("%d", 1:p), rho=0.5, g=NULL, tol=0.001, verbose=FALSE) {
  p <- as.integer(p)
  d <- as.integer(d)

  if (!is.null(g)) {
    if (nrow(as(g, "matrix")) != p)
      stop("'p' does not match the number of vertices in 'g'\n")
  }

  if (rho <= -1/(p-1) || rho >= 1)
    stop("'rho' should be a real number such that -1/(p-1) < 'rho' < 1.")

  sim <- list()
  for (i in 1:n) {
    sim.g <- g
    if (is.null(sim.g))
      sim.g <- qpRndGraph(p=p, d=d, labels=labels, return.type="graphBAM", verbose=verbose)

    sim.sigma <- qpG2Sigma(sim.g, rho=rho, tol=tol, verbose=verbose)

    sim[[i]] <- UGgmm(g=sim.g, mean=rep(0, p), sigma=sim.sigma)
  }

  if (n == 1)
    sim <- sim[[1]]

  sim
}



## function: rHMgmm
## purpose: random generation of an homogeneous mixed graphical Markov model (HMGMM)
## parameters: n - number of HMGMMs to simulate, by default only one
##             pI - number of discrete r.v.'s
##             pY - number of continuous r.v.'s
##             d - degree of every vertex in the undirected graph
##             labels - vertex labels corresponding to r.v. names
##             rho - mean marginal correlations for the present pure continuous edges
##             a - mean difference in the continuous r.v.'s forming part of a mixed edge
##             dLevels - number of discrete levels
##             g - given graph, if only the covariance matrix needs to be simulated, by
##                 default is NULL which implies that a random graph is also simulated
##             tol - maximum tolerance used by the matrix completion algorithm when
##                   simulating a random covariance matrix
##             verbose - show progress on the calculations, specially useful for
##                       high dimensional UGGMs with p > 500
## return: an HMgmm object containing the simulated HMGMM if n=1 and a list of HMgmm objects
##         with the simulated HMGMMs if n > 1

oldrHMgmm <- function(n=1L, pI=1L, pY=4L, d=2L, labels=c(sprintf("I%d", 1:pI), sprintf("Y%d", 1:pY)),
                      rho=0.5, a=1, dLevels=2L, g=NULL, tol=0.001, verbose=FALSE) {
  if (pI < 1)
    stop("pI should be equal or larger than 1.")
  if (pY < 1)
    stop("pY should be equal or larger than 1.")

  p <- as.integer(pI + pY)

  if (!is.null(labels))
    if (length(labels) != p)
      stop("'labels' should be either NULL or 'pI' + 'pY' variable names.")

  if (length(a) == 1)
    a <- rep(a, pY)
  else if (length(a) != pY)
    stop("the vector of additive effects 'a' should contain 'pY' values.")

  if (!is.null(g)) {
    if (class(g) == "graphBAM") {
      if (is.na(match("type", names(nodeData(g, nodes(g))))))
        stop("When 'g' is a graphBAM object, it should contain node type (discrete, continuous) attribute (see graph::nodeData)\n")
    } else if (class(g) == "matrix" || length(grep("Matrix", class(g))) > 0) {
      dimg <- dim(g)
      if (dimg[1] != dimg[2])
        stop("'g' is not an squared matrix.")

      if (dimg[1] != p)
        stop("'g' should be a p x p squared matrix with p = #discrete variables  +  #continuous variables.")

      if (!isSymmetric(g))
        stop("'g' is not symmetric.")

      if (class(g[1, 1]) == "integer" || class(g[1, 1]) == "numeric") {
        if (verbose)
          warning("coercing input numeric adjacency matrix 'g' to a logical adjacency matrix\n")

        g <- g != 0
      }
      g <- as(g, "matrix") ## coerce to a classical Matrix since subsetting is still not possible with Matrix-matrices
      if (!is.null(colnames(g))) ## column names of input adjacency matrix override 'labels' argument
        labels <- colnames(g)

      if (is.null(labels))
        labels <- c(sprintf("I%02d", 1:pI), sprintf("Y%02d", 1:pY))

      if (is.null(colnames(g)))
        rownames(g) <- colnames(g) <- labels

      df <- data.frame(from=labels[row(g)[upper.tri(g) & g]],
                       to=labels[col(g)[upper.tri(g) & g]],
                       weight=rep(1, sum(g)/2))
      g <- graphBAM(df, edgemode = "undirected", nodes=labels)
      nodeDataDefaults(g, "type") <- "continuous"
      nodeData(g, labels[1:pI], "type") <- "discrete"
      ## the following two lines are necessary since graphBAM() reorders vertices alphabetically
      labels <- nodes(g)
      vtype <- factor(unlist(nodeData(g, nodes(g), "type")))

      if (is.null(names(a)))
        names(a) <- labels[vtype == "continuous"]
      else {
        if (any(is.na(match(names(a), labels[vtype == "continuous"]))))
          stop("some continuous variable names in 'a' are not part of the variable names in 'g'.")
      }
    } else
      stop("'g' should be either a graphBAM object or an adjacency matrix.\n")

    vtype <- factor(unlist(nodeData(g, nodes(g), "type")))
    ed <- edges(g)[vtype == "discrete"]
    if (any(sapply(ed, function(x, vt) any(vt[x] == "discrete"), vtype)))
      stop("'g' cannot contain edges between vertices of discrete variables.")
  }

  if (rho <= -1/(pY-1) || rho >= 1)
    stop("'rho' should be a real number such that -1/(pY-1) < 'rho' < 1.")

  if (any(dLevels > 2))
    stop("Only binary variables can be used at the moment.")

  if (any(dLevels < 1))
    stop("Discrete variables should have at least two levels specified in argument 'dLevels'.")

  sim <- list()
  for (i in 1:n) {
    sim.g <- g
    if (is.null(sim.g)) {
      sim.g <- qpRndGraph(p=p, d=d, labels=labels, exclude=1:pI, return.type="graphBAM", verbose=verbose)
      nodeDataDefaults(sim.g, "type") <- "continuous"
      nodeData(sim.g, labels[1:pI], "type") <- "discrete"
    }

    vtype <- factor(unlist(nodeData(sim.g, nodes(sim.g), "type"), use.names=FALSE))
    Y <- graph::nodes(sim.g)[vtype == "continuous"]
    I <- graph::nodes(sim.g)[vtype == "discrete"]

    if (is.null(names(a)))
      names(a) <- Y

    sim.sigma <- qpG2Sigma(g=subGraph(Y, sim.g), rho=rho, verbose=verbose)

    sim[[i]] <- HMgmm(g=sim.g, dLevels=dLevels, a=a, rho=rho, sigma=sim.sigma)
  }

  if (n == 1)
    sim <- sim[[1]]

  sim
} 


## function: qpRndHMGM
## purpose: builds a random homogeneous mixed graphical Markov model
##          and for every vertex its boundary <= d
## parameters: nDiscrete - number of discrete variables
##             nContinuous - number of continuous variables
##             d - degree of every vertex
##             mixedIntStrength - strength of the mixed interactions
##             rho - marginal correlation of the quadratic interactions
##             G - input graph, if given
## return: a list with the HMGM

qpRndHMGM <- function(I=1, Y=3, d=2, rho=0.5, mudif=NULL, eta2=0.1, Ilevels=2, g=NULL, verbose=TRUE) {
  if (length(I) == 1)
    I <- paste0("I", 1:as.integer(I))

  if (length(Y) == 1)
    Y <- paste0("Y", 1:as.integer(Y))

  pI <- length(I)
  pY <- length(Y)
  p <- pI + pY

  if (pI == 0)
    stop("At least one discrete variable should be specified in argument I.")

  if (pY == 0)
    stop("At least one continuous variable should be specified in argument Y.")

  if (length(Ilevels) == 1)
    Ilevels <- rep(Ilevels, pI)
  else {
    if (length(Ilevels) != pI)
      stop("Specified discrete levels in argument Ilevels do not match the number of discrete variables in argument I.")
  }

  if (any(Ilevels > 2))
    stop("This function cannot handle yet discrete variables with more than two levels.")

  if (any(Ilevels < 2))
    stop("Discrete variables should have at least two levels specified in argument Ilevels.")

  if (!is.null(mudif)) {
    eta2 <- NULL
    if (length(mudif) == 1)
      mudif <- rep(mudif, pY)
    else if (length(mudif) != pY)
      stop("Specified mean differences in argument mudif do not match the number of continuous variables in argument Y.")
  }

  if (!is.null(eta2)) {
    if (length(eta2) == 1)
      eta2 <- rep(eta2, pY)
    else if (length(eta2) != pY)
      stop("Specified standardized effect sizesin argument eta2 do not match the number of continuous variables in argument Y.")
  }

  if (is.null(g)) {
    if (verbose)
      message(sprintf("Generating a random d-regular graph G of constant degree %d.", d))
    g <- qpRndGraph(p=p, d=d, exclude=1:pI, verbose=verbose)
    rownames(g) <- colnames(g) <- c(I, Y)
  } else {
    dimg <- dim(g)
    if (dimg[1] != dimg[2])
      stop("g is not an squared matrix")

    if (dimg[1] != p)
      stop("g should be a p x p squared matrix with p = #discrete variables I + #continuous variables Y")

    if (is.null(dimnames(g)))
      dimnames(g) <- list(c(I, Y), c(I, Y))

    if (!identical(rownames(g), c(I, Y)) || !identical(colnames(g), c(I, Y)))
      stop("row and column names of g should be identical to the discrete and continuous variable names in arguments I and Y")
  }

  ## if (any(rowSums(as.matrix(g[I, Y])) > 1)) ## rowSums() does not seem to work with Matrix matrices
  if (any(Matrix::rowSums(g[I, Y]) > 1)) ## rowSums() does not seem to work with Matrix matrices
    stop("Procedure not yet prepared for connecting a discrete variable to more than one continuous variable.")

  ## if (any(rowSums(as.matrix(g[Y, I])) > 1)) ## rowSums() does not seem to work with Matrix matrices
  if (any(Matrix::rowSums(g[Y, I]) > 1)) ## rowSums() does not seem to work with Matrix matrices
    stop("Procedure not yet prepared for connecting a continuous variable to more than one discrete variable.")

  if (verbose)
    message("Generating a random covariance matrix whose inverse has a pattern of zeroes defined by the missing edges in G.")

  ## simulate the covariance matrix Sigma
  Sigma <- qpG2Sigma(g[Y, Y], rho=rho, verbose=verbose)
  rownames(Sigma) <- colnames(Sigma) <- Y

  if (is.null(mudif)) { ## simulate mean differences according to specified standardized effect sizes eta2
    m <- Sigma
    m <- t(sapply(1:pY, function(i, x, v, h) -h[i]*(x[i, ]/v)^2,
                  m, diag(Sigma), eta2))
    diag(m) <- 1 - eta2
    mudif2 <- solve(m)%*%(eta2*diag(Sigma))
    mudif <- 2*as.vector(sqrt(mudif2))
    names(mudif) <- Y
  } else ## use specified mean differences and calculate resulting standardized effect sizes eta2
    eta2 <- sapply(1:pY, function(i) ((mudif[i]/2)^2) /
                       (Sigma[i, i] + sum(sapply(1:pY, function(j) ((Sigma[i,j] / Sigma[j,j])^2) * (mudif[j]/2)^2))))
  
  list(I=I, Y=Y, Ilevels=Ilevels, g=g, d=d, Sigma=Sigma, rho=rho, mudif=mudif, eta2=eta2)
}



## function: qpSampleFromHMGM
## purpose: samples synthetic data from a homogeneous mixed graphical Markov model
## parameters: n - number of observations
##             hmgm - model as generated by the function qpRndHMGM()
## return: the sampled synthetic data

qpSampleFromHMGM <- function(n=10, hmgm=qpRndHMGM()) {

  Ilevels <- hmgm$Ilevels
  pI <- length(hmgm$I)
  pY <- length(hmgm$Y)
  p <- pI + pY

  X <- matrix(1, nrow=n, ncol=p, dimnames=list(1:n, c(hmgm$I, hmgm$Y)))

  ## simulate discrete data uniformly at random
  if (length(unique(Ilevels)) == 1)
    X[, 1:pI] <- sample(1:Ilevels[1], size=n*pI, replace=TRUE)
  else
    stop("This function cannot simulate yet data from discrete variables with different number of levels.")

  ## we use the same data matrix X to store the intermediate values required
  ## during the simulation process.
 
  ## we start by setting on the continuous variables what data points are
  ## associated to a particular discrete level, i.e., from what Gaussian
  ## distribution a particular data points must be sampled from
  whYxI <- which(rowSums(as.matrix(hmgm$g[hmgm$Y, hmgm$I])) > 0) ## which Yi are connected to I rv's
  X[, whYxI+pI] <- sapply(whYxI, function(Yi, X, pI, mod) {
                            whIxYi <- which(mod$g[Yi+pI, mod$I])
                            apply(X[, whIxYi, drop=FALSE], 1,
                                 function(i) sum((i-1)*2^((length(i)-1):0)))
                          }, X, pI, hmgm)

  ## calculate canonical parameter h(i) from the homogeneous mixed graphical model
  ## by now this only works with at most one discrete variable associated to a continuous one
  whIxY <- which(rowSums(as.matrix(hmgm$g[hmgm$I, hmgm$Y, drop=FALSE])) > 0) ## which Ii are connected to Y rv's
  for (Ii in whIxY) {
    whYxIi <- which(hmgm$g[Ii, hmgm$Y])
    if (length(whYxIi) == 1)
      X[, pI+whYxIi] <- X[, pI+whYxIi] * hmgm$mudif[whYxIi] / hmgm$Sigma[whYxIi, whYxIi]
    else {
      difh <- as.matrix(solve(hmgm$Sigma[whYxIi, whYxIi])) %*% hmgm$mudif[whYxIi]
      for (j in whYxIi)
        X[, j+pI] <- X[, j+pI] * difh[j, ]
    }
  }

  ## mu = Sigma x h where corresponding h(i) are stored in X[, Y] for convenience
  X[, hmgm$Y] <- t(as.matrix(hmgm$Sigma) %*% t(X[, hmgm$Y]))

  ## sample observations from the homogenous mixed graphical model
  xtab <- tapply(1:n, apply(X[, whIxY, drop=FALSE], 1, function(i) paste(i, collapse="")))
  xtab <- split(as.data.frame(X[, hmgm$I]), xtab)
  for (i in 1:length(xtab)) {
    li <- xtab[[i]]
    which_n <- as.numeric(rownames(li))
    m <- X[which_n, hmgm$Y, drop=FALSE]
    X[which_n, hmgm$Y] <- mvtnorm::rmvnorm(length(which_n), mean=m[1, ],
                                           sigma=as.matrix(hmgm$Sigma))
  }

  return(X)
}
