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



## function: rUGgmm
## purpose: random generation of an undirected Gaussian graphical Markov model (UGGMM)
## parameters: n - number of UGGMMs to simulate, by default only one
##             p - number of r.v.'s
##             d - degree of every vertex in the undirected graph
##             rho - mean marginal correlations for the present edges
##             g - given graph, if only the covariance matrix needs to be simulated, by
##                 default is NULL which implies that a random graph is also simulated
##             tol - maximum tolerance used by the matrix completion algorithm when
##                   simulating a random covariance matrix
##             verbose - show progress on the calculations, specially useful for
##                       high dimensional UGGMs with p > 500
## return: an UGgmm object containing the simulated UGGMM if n=1 and a list of UGgmm objects
##         with the simulated UGGMMs if n > 1

rUGgmm <- function(n=1L, p=5L, d=2L, labels=sprintf("%02d", 1:p), rho=0.5, g=NULL, tol=0.001, verbose=FALSE) {
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

rHMgmm <- function(n=1L, pI=1L, pY=4L, d=2L, labels=c(sprintf("I%02d", 1:pI), sprintf("Y%02d", 1:pY)),
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

## genes: # of genes or vector of cM positions, each element corresponding to distinct gene
## cis: fraction of genes with eQTL in cis
## trans: vector of integer numbers, each element corresponding to a trans-eQTL and its value
##        is the number of genes linked to this eQTL
## cisr: cis radius, maximum distance in cM that defines a QTL in cis wrt to a gene
## d2m: distance of every gene or eQTL to a marker, default 0 implies eQTLs and genes are located at markers,
##      and therefore, the maximum number of eQTL is bounded by the number of markers
reQTLcross <- function(n=1L, map, type=c("bc"), genes=100, cis=1, trans=c(), cisr=1, d2m=0,
                       d=2, g=NULL, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
  type <- match.arg(type)

  if (class(n) != "numeric" && class(n) != "integer")
    stop("argument 'n' should contain the number of eQTLcross models to simulate (default=1).")

  if (cis < 0 || cis > 1)
    stop("argument 'cis' should be a real number between 0 and 1.")

  nm <- qtl::totmar(map)                 ## total number of markers
  nmbychr <- qtl::nmar(map)              ## number of markers per chromosome
  csnmbychr <- cumsum(nmbychr)           ## cumulative sum of the number of markers per chromosome
  mcmloc <- unlist(map, use.names=FALSE) ## location of markers in cM
  cmlenbychr <- sapply(map, max)         ## chromosome length in cM
  cscmlenbychr <- cumsum(cmlenbychr)     ## cumulative sum of chromosome length in cM

  pY <- length(genes)
  if (length(genes) == 1)
    pY <- genes

  if (pY > nm && d2m == 0)
    stop("more genes than markers. Either, increase marker density in the genetic map, increase d2m or decrease the number of genes.")

  Y <- sprintf("g%d", 1:pY)

  chr.genes <- NA
  ## simulate gene locations in cM
  if (length(genes) == 1) {
    nocis <- FALSE
    i <- 0
    ## enforce genes being located at least 2 x cisr cM apart
    while (!nocis && i < 10) {
      tmpg <- sample(1:nm, size=genes, replace=FALSE)
      chr.genes <- sapply(tmpg, function(i, cs) sum(cs < i)+1, csnmbychr)
      tmpg <- mcmloc[tmpg] + d2m
      nocis <- sapply(split(tmpg, chr.genes), function(gxc, cisr) {
                                             nocis <- TRUE
                                             if (length(gxc) > 1)
                                               nocis <- all(combn(gxc, 2, function(x) abs(x[1]-x[2])) > 2*cisr)
                                             nocis
                                           }, cisr)
      nocis <- all(nocis)
      if (nocis)
        genes <- tmpg

      i <- i + 1
    }

    if (!nocis)
      stop("impossible to simulate genes. Either decrease cisr, decrease the number of genes, or increase marker density in the genetic map.")
  }

  ## build gene annotation matrix
  n.genes <- length(genes)
  genes <- cbind(chr.genes, genes)
  genes <- genes[order(genes[, 1], genes[, 2]), ]
  colnames(genes) <- c("chr", "location")
  rownames(genes) <- Y

  n.cisQTL <- floor(n.genes * cis)       ## number of cisQTL
  n.transQTL <- length(trans)            ## number of transQTL

  if (n.cisQTL+n.transQTL > nm && d2m == 0)
    stop("more eQTL than markers. Either, increase marker density in the genetic map, increase d2m or decrease the number of eQTL.")

  if ((class(a) == "numeric" || class(a) == "integer") && length(a) > 1 && length(a) != n.cisQTL + n.transQTL)
    stop(sprintf("argument 'a' contains %d values of eQTL additive effects while arguments 'genes', 'cis' and 'trans' determine a total number of %d eQTL.", length(a), n.cisQTL+n.transQTL))

  if (class(a) == "function" && length(formals(a)) != 1)
    stop("when argument 'a' is a function it should contain one argument taking the number of eQTL.")

  ## function to search markers (m) in cis to a gene (g) within a radius (r)
  cism <- function(markers, gene, radius) which(markers >= gene-radius & markers <= gene+radius)

  pI <- n.cisQTL + n.transQTL
  I <- sprintf("QTL%d", 1:pI)

  dLevels <- switch(type, bc=2, NA)

  sim <- list()
  for (i in 1:n) {
    ## simulate cis-QTL associations
    cisQTL <- matrix(NA, nrow=0, ncol=3)
    cisQTLgenes <- sample(1:n.genes, size=n.cisQTL, replace=FALSE) ## which genes should this cis-markers associated to?
    cisQTLgenes <- split(cisQTLgenes, names(map)[genes[cisQTLgenes, "chr"]])
    for (chr in names(cisQTLgenes)) {
      loc.genes <- genes[cisQTLgenes[[chr]], "location"]
      markers <- map[[chr]] + d2m
      allcm <- sapply(loc.genes, function(gene, markers, cisr) cism(markers, gene, cisr),
                      markers, cisr, simplify=FALSE) ## all cis-markers
      cm <- c(1, 1)
      j <- 1
      while (any(duplicated(cm)) && j < 10) {
        cm <- sapply(allcm, function(x) { if (length(x) > 1) x <- sample(x, size=1) ; x}) ## select one cis-marker
        j <- j + 1
      }
      if (any(duplicated(cm)))
        stop("impossible to simulate cis-eQTL. Either decrease cisr or increase marker density in the genetic map.")

      cisQTL <- rbind(cisQTL, cbind(rep(match(chr, names(map)), times=length(cm)), markers[cm], cisQTLgenes[[chr]]))
    }

    ## simulate trans-QTL associations
    transQTL <- matrix(NA, nrow=0, ncol=3)
    transgenes <- setdiff(1:n.genes, cisQTL[, 3])
    if (sum(trans) > length(transgenes))
      stop("not enough genes to simulate trans-QTL associations. Either decrease the number of trans-QTL associations or decrease 'cis'.")

    for (ng in trans) {
      tmpmap <- map
      transQTLgenes <- sample(transgenes, size=ng, replace=FALSE)
      transgenes <- setdiff(transgenes, transQTLgenes) ## removed sampled trans-genes
      transQTLgenes <- split(transQTLgenes, names(tmpmap)[genes[transQTLgenes, "chr"]])
      for (chr in names(transQTLgenes)) {
        loc.genes <- genes[transQTLgenes[[chr]], "location"]
        allcm <- sapply(loc.genes, function(g, m, cisr) cism(m, g, cisr), tmpmap[[chr]]+d2m, cisr) ## all cis-markers
        tmpmap[[chr]] <- tmpmap[[chr]][-allcm] ## remove all cis-markers
      }
      tmpmcmloc <- unlist(tmpmap, use.names=FALSE) + d2m
      tmpnmbychr <- qtl::nmar(tmpmap)
      tmpcsnmbychr <- cumsum(tmpnmbychr)
      tm <- sample(sum(tmpnmbychr), size=1, replace=FALSE)
      chr.tm <- sum(tmpcsnmbychr < tm) + 1
      transQTL <- rbind(transQTL, cbind(rep(chr.tm, times=ng), rep(tmpmcmloc[tm], times=ng),
                                        unlist(transQTLgenes, use.names=FALSE)))
    }
    
    ## simulate gene network
    sim.g <- g
    if (is.null(sim.g))
      sim.g <- qpRndGraph(p=pY, d=d, labels=Y, return.type="graphBAM", verbose=verbose)
    edges <- graph::edges(sim.g)
    edges <- cbind(rep(names(edges), times=sapply(edges, length)),
                   unlist(edges, use.names=FALSE))
    edges <- unique(t(apply(edges, 1, sort)))
    edges <- cbind(match(edges[, 1], Y), match(edges[, 2], Y))

    ## simulate conditional covariance matrix
    sim.sigma <- qpG2Sigma(g=sim.g, rho=rho, tol=tol, verbose=verbose)

    ## simulate additive effect in QTL
    qtl <- rbind(cisQTL, transQTL)

    if ((class(a) == "numeric" || class(a) == "integer") && length(a) == 1)
      a <- rep(a, times=nrow(qtl))
    else if (class(a) == "function")
      a <- a(nrow(qtl))

    qtl <- cbind(qtl, a)

    sim[[i]] <- eQTLcross(map, genes=genes, model=qtl, type=type, nGenes=pY,
                          geneNetwork=edges, rho=rho, sigma=sim.sigma)
  }

  if (n == 1)
    sim <- sim[[1]]

  sim
}

## overload sim.cross() from the qtl package to enable simulateing eQTL data from
## an experimental cross

sim.cross <- function(map, model, ...) UseMethod("sim.cross", model)
sim.cross.default <- function(map, model, ...) qtl::sim.cross(map, model, ...)
sim.cross.matrix <- function(map, model, ...) qtl::sim.cross(map, model, ...)
setMethod("sim.cross", c(map="map", model="matrix"), sim.cross.matrix)

sim.cross.eQTLcross <- function(map, model, n.ind=100, ...) {
            crossModel <- unique(alleQTL(model)[, c("chrom", "location")])
            crossModel <- cbind(crossModel, dummya=rep(1, nrow(crossModel)))
            cross <- qtl::sim.cross(map=map, model=crossModel, type=model$type, n.ind=n.ind, ...)

            ## we use the same data matrix X to store the mean values employed
            ## during the simulation process.
            I <- model@model$I
            Y <- model@model$Y
            stopifnot(identical(I, colnames(cross$qtlgeno)))
            cross$pheno <- qpgraph:::calculateCondMean(model@model, cross$qtlgeno) 
            rownames(cross$pheno) <- 1:n.ind
            colnames(cross$pheno) <- Y
            cross$pheno <- as.data.frame(cross$pheno)

            YxI <- Y[which(sapply(graph::edges(model@model$g)[Y], function(xYi, vt) sum(vt[xYi] == "discrete"), model@model@vtype) > 0)]
            xtab <- tapply(1:n.ind, apply(cross$pheno[, YxI, drop=FALSE], 1, function(i) paste(i, collapse="")))
            xtab <- split(as.data.frame(cross$qtlgeno[, I]), xtab)
            for (i in 1:length(xtab)) {
              li <- xtab[[i]]
              which_n <- as.numeric(rownames(li))
              cross$pheno[which_n, Y] <- mvtnorm::rmvnorm(length(which_n), mean=as.numeric(cross$pheno[which_n[1], Y]),
                                                          sigma=as.matrix(model@model$sigma))
            }

            cross
          }
setMethod("sim.cross", signature(map="map", model="eQTLcross"), sim.cross.eQTLcross)

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

qpSimulateEqtl <- function(nContinuous=10, d=2, nHotspots=0, h2=0.5, rho=0.5, G=NULL, verbose=TRUE) {

  Delta <- paste("D", 1:nContinuous, sep="")
  Gamma <- paste("C", 1:nContinuous, sep="")
  if (is.null(G)) {
    G <- qpRndGraph(p=length(Gamma), d=d)
    rownames(G) <- colnames(G) <- Gamma
  } else {
    Gamma <- colnames(G)
  }
  
  nDiscreteLevels <- 2
  eQTLs <- as.list(1:nContinuous) ## stores for each gene to what discrete var is linked
  names(eQTLs) <- Gamma
  if (nHotspots > 0) {
    hotspots <- sample(1:nContinuous, size=nHotspots, replace=FALSE)
    x <- runif(nContinuous, min=0, max=1)
    for (i in 1:nContinuous) {
      if (length(intersect(i, hotspots)) == 0) {
	      if (x[i] < 0.5) {
	        whotspot <- ifelse(nHotspots > 1, sample(hotspots, size=1), hotspots)
	        eQTLs[[i]] <- whotspot
	      }
      }
    }
  }

  ## Sigma
  Sigma <- qpG2Sigma(G, rho=rho, verbose=verbose)
  rownames(Sigma) <- colnames(Sigma) <- Gamma

  ## additive effect as function of heritability (h2)  
  m <- Sigma  
  m <- t(apply(m, 1, function(x, v, h) -h*(x/v)^2, diag(as.matrix(Sigma)), h2))
  diag(m) <- 1 - h2
  
  a2 <- solve(m)%*%(h2*diag(as.matrix(Sigma)))
  a <- 2*as.vector(sqrt(a2))
  names(a) <- Gamma
 
  list(Delta=Delta, Gamma=Gamma, G=G, nDiscreteLevels=nDiscreteLevels, eQTLs=eQTLs, Sigma=Sigma, h2=h2, a=a)
}


qpSampleFromEqtl <- function(n=10, hmgm=qpSimulateEqtl()) {
  require(mvtnorm)

  nDiscreteLevels <- hmgm$nDiscreteLevels
  nDiscrete <- length(hmgm$Delta)
  nContinuous <- length(hmgm$Gamma)

  sampleData <- matrix(1, nrow = n, ncol = (nDiscrete + nContinuous), dimnames = list((1:n), c(hmgm$Delta, hmgm$Gamma)))

  ## simulate discrete data uniformly random  
  sampleData[, 1:nDiscrete] <- sample(1:nDiscreteLevels, size=n*nDiscrete, replace=TRUE)  
  
  ## mark associated levels in continuous observations
  for (i in 1:n) {
    for (j in 1:length(hmgm$eQTLs)) {
      g <- names(hmgm$eQTLs)[j]
      Delta <- hmgm$eQTLs[[j]] ## under no epistasis this will be always a singleton
      sampleData[i, g] <- sum((sampleData[i, Delta] - 1) * 2^((length(Delta)-1):0))
    }
  }

  ## calculate h in the mixed graphical model as the additive effect
  for (i in 1:nContinuous) { ## this is in fact going through each discrete variable (marker)
    GammabyDelta <- names(which(unlist(lapply(hmgm$eQTLs, function(x, i) {x == i}, i))))
    if (length(GammabyDelta) > 0) {
      if (length(GammabyDelta) == 1) {
	      sampleData[ ,GammabyDelta] <- sampleData[ ,GammabyDelta]*hmgm$mudif[GammabyDelta]/hmgm$Sigma[GammabyDelta, GammabyDelta]
      } else {
        difh <- solve(hmgm$Sigma[GammabyDelta, GammabyDelta])%*%hmgm$mudif[GammabyDelta]
        for (j in GammabyDelta)
	        sampleData[ ,j] <- sampleData[ ,j]*difh[j, ]          
      }
    }
  }   
 
    
  ## multiply Sigma * h to obtain mu  
  sampleData[, hmgm$Gamma] <- t(as.matrix(hmgm$Sigma) %*% t(sampleData[, hmgm$Gamma]))
  
  ## xtab <- tapply(1:n, as.data.frame(sampleData[, hmgm$Delta]))

  xtab <- tapply(1:n, apply(sampleData[, hmgm$Delta], 1, function(x) paste(x, collapse="")))
  xtab <- split(as.data.frame(sampleData[, hmgm$Delta]), xtab)
  for (i in 1:length(xtab)) {
    li <- xtab[[i]]
    which_n <- as.numeric(rownames(li))    
    m <- sampleData[which_n, hmgm$Gamma, drop=FALSE]
    sampleData[which_n, hmgm$Gamma] <- rmvnorm(length(which_n), mean=m[1, ], 
                                               sigma=as.matrix(hmgm$Sigma))
  }    
  
  return(sampleData)
}

