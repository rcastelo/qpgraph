setValidity("HMgmm",
            function(object) {
              valid <- TRUE

              if (graph::edgemode(object@g) != "undirected")
                valid <- "the underlying graph 'g' should have 'edgemode=\"undirected\" even though 'g' is treated as semi-directed'."

              if (class(valid) == "logical" && length(object@a) != sum(object@vtype == "continuous"))
                valid <- "the vector of additive effects 'a' should contain as many values as continuous variables."

              if (class(valid) == "logical" && graph::numNodes(object@g) != object$pI+object$pY)
                valid <- "the number of vertices in 'g' does not match the number of variables."

              p <- (d <- dim(object$sigma))[1]
              if (class(valid) == "logical" && p != d[2])
                valid <- "'sigma' should be a squared matrix."

              if (class(valid) == "logical" && p != sum(object@vtype == "continuous"))
                valid <- "the dimension of sigma should match the number of continuous variables."

              if (class(valid) == "logical" && !isSymmetric(as(object$sigma, "matrix")))
                valid <- "'sigma' should be symmetric."

              if (class(valid) == "logical" && !identical(colnames(object$sigma), object$Y))
                valid <- "column names in 'sigma' should be identical to continuous variable names in 'Y'."

              if (class(valid) == "logical" && !identical(unlist(nodeData(object@g, nodes(object@g), "type"),
                                                                 use.names=FALSE),
                                                          as.character(object@vtype)))
                valid <- "the type of vertices in 'vtype' does not match the node type information in 'g'."

              if (class(valid) == "logical" && !identical(nodes(object@g), names(object@vtype)))
                valid <- "'vtype' should be a named vector whose names match the vertex labels in 'g'."

              if (class(valid) == "logical" && !identical(nodes(object@g)[object@vtype == "continuous"], names(object@a)))
                valid <- "'a' should be a named vector whose names match the vertex labels in 'g'."

              valid
            })


## constructor methods
setMethod("HMgmm", signature(g="missing"),
          function(g, pI=1, pY=4, dLevels=2L) {
            dVertexLabels <- sprintf("I%d", seq(along=rep(1, pI))) ## to deal with pI=0
            cVertexLabels <- sprintf("Y%d", seq(along=rep(1, pY))) ## to deal with pY=0
            g <- graph::graphBAM(as.data.frame(matrix(NA, nrow=0, ncol=3,
                                                      dimnames=list(NULL, c("from", "to", "weight")))),
                                 nodes=c(dVertexLabels, cVertexLabels))
            nodeDataDefaults(g, "type") <- "continuous"
            nodeData(g, dVertexLabels, "type") <- "discrete"
            vtype <- factor(unlist(nodeData(g, nodes(g), "type"), use.names=FALSE))

            a <- rep(NA, pY)
            names(a) <- nodes(g)[vtype == "continuous"]

            sigma <- new("dspMatrix", Dim=as.integer(c(pY, pY)),
                        Dimnames=list(cVertexLabels, cVertexLabels),
                        x=diag(pY)[upper.tri(diag(pY), diag=TRUE)])

            HMgmm(g=g, dLevels=dLevels, a=a, rho=0, sigma=sigma)
          })

setMethod("HMgmm", signature(g="graphBAM"),
          function(g, dLevels=2L, a=0, rho=0.5,
                   sigma=diag(sum(unlist(graph::nodeData(g, graph::nodes(g), "type"), use.names=FALSE) == "continuous"))) {

            vtype <- factor(unlist(graph::nodeData(g, graph::nodes(g), "type"), use.names=FALSE))
            names(vtype) <- nodes(g)

            pI <- sum(vtype == "discrete")
            pY <- sum(vtype == "continuous")
            I <- nodes(g)[vtype == "discrete"]
            Y <- nodes(g)[vtype == "continuous"]

            maskYxI <- sapply(graph::nodes(g)[vtype == "continuous"],
                              function(v, e, vt) any(vt[e[[v]]] == "discrete"),
                              graph::edges(g), vtype)

            if (length(a) == 1) {
              tmp <- rep(NA, pY)
              names(tmp) <- Y
              tmp[maskYxI] <- rep(a, sum(maskYxI))
              a <- tmp
            }

            if (is.null(names(a)))
              stop("when 'a' is a vector it should be named vector of linear additive effects with continuous variable names.")
            else if (any(is.na(match(names(a), Y))))
              stop("some variable names in 'a' do not form part of the continuous variables in the model.")
            else if (length(a) < pY) {
              tmp <- rep(NA, pY)
              names(tmp) <- Y
              tmp[names(a)] <- a
              a <- tmp
            }

            new("HMgmm", pI=pI, pY=pY, g=g, vtype=vtype, dLevels=dLevels, a=a,
                rho=rho, sigma=as(sigma, "dspMatrix"), eta2=new.env(parent=emptyenv()))
          })

## constructor simulation methods
setMethod("rHMgmm", signature(n="markedGraphParam", g="missing"),
          function(n, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g=n, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="missing", g="markedGraphParam"),
          function(n, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="numeric", g="markedGraphParam"),
          function(n=1, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=as.integer(n), g, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="integer", g="markedGraphParam"),
          function(n=1L, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            pI <- g@pI
            pY <- g@pY

            if (pI < 1)
              stop("pI should be equal or larger than 1.")
            if (pY < 1)
              stop("pY should be equal or larger than 1.")

            p <- as.integer(pI + pY)

            labels <- c(g@Ilabels, g@Ylabels)

            if (!is.null(labels))
              if (length(labels) != p)
                stop("when set, vertex labels should contain 'pI' + 'pY' variable names.")

            if (rho <= -1/(pY-1) || rho >= 1)
              stop("'rho' should be a real number such that -1/(pY-1) < 'rho' < 1.")

            if (any(dLevels > 2))
              stop("Only binary variables can be used at the moment.")

            if (any(dLevels < 1))
              stop("Discrete variables should have at least two levels specified in argument 'dLevels'.")

            if (is(g, "dRegularMarkedGraphParam")) {
              if (g@pI > 1)
                g@exclude <- 1:g@pI ## assuming we put consider always discrete vertices before continuous ones
            }

            sim <- list()
            for (i in 1:n) {
              sim.g <- rgraphBAM(g)

              if (is(g, "erMarkedGraphParam")) { ## remove edges between discrete vertices with the ER model
                dEdges <- edges(sim.g)[g@Ilabels]
                ndEdges <- sapply(dEdges, function(x, I) sum(x %in% I), g@Ilabels)
                if (any(ndEdges > 0)) {
                  dEdges <- dEdges[ndEdges > 0]
                  nEdges <- nEdges[nEdges > 0]
                  dEdges <- lapply(dEdges, function(x, I) x[x %in% I], g@Ilabels)
                  dEdges <- cbind(rep(names(dEdges), nEdges), unlist(dEdges, use.names=FALSE))
                  sim.g <- graph::removeEdge(from=dEdges[, 1], to=dEdges[, 2], sim.g)
                }
              }
              nodeDataDefaults(sim.g, "type") <- "continuous"
              nodeData(sim.g, g@Ilabels, "type") <- "discrete"

              sim.sigma <- qpG2Sigma(g=subGraph(g@Ylabels, sim.g), rho=rho, verbose=verbose)
              sim.sigma <- sim.sigma[g@Ylabels, g@Ylabels] ## put back rows and columns into the original variable order
                                                           ## since 'subgraph()' re-orders nodes alphabetically (sigh!)

              sim[[i]] <- HMgmm(g=sim.g, dLevels=dLevels, a=a, rho=rho, sigma=sim.sigma)
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

setMethod("rHMgmm", signature(n="matrix", g="missing"),
          function(n, g, I=1, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g=n, I, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="missing", g="matrix"),
          function(n, g, I=1, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g, I, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="numeric", g="matrix"),
          function(n=1, g, I=1, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=as.integer(n), g, I, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="integer", g="matrix"),
          function(n=1, g, I=1, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            p <- (d <- dim(g))[1]
            if (p != d[2] && d[2] != 2)
              stop("If 'g' is a matrix it should be either a squared symmetric adjacency matrix or a two-column matrix with vertex pairs in the rows defining the edge set of a marked graph.")

            df <- as.data.frame(matrix(NA, nrow=0, ncol=3, dimnames=list(NULL, c("from", "to", "weight"))),
                                                                stringsAsFactors=FALSE)
            vlabels <- Ilabels <- Ylabels <- NULL
            if (p == d[2]) {
              if (!isSymmetric(g))
                stop("'g' is not a symmetric matrix\n")

              if (class(g[1, 1]) == "integer" || class(g[1, 1]) == "numeric") {
                if (verbose)
                  warning("coercing input numeric adjacency matrix 'g' to a logical adjacency matrix\n")

                g <- g != 0
              }

              vlabels <- colnames(g)
              if (is.null(vlabels)) {
                if (is.character(I))
                  stop("if 'I' specifies character vertex labels, then 'g' should contain row and column names of all of the vertex labels.")
                Ilabels <- sprintf("I%d", seq(along=I))
                Ylabels <- sprintf("Y%d", seq(along=setdiff(1:p, I)))
                vlabels <- sprintf("X%d", 1:p)
                vlabels[setdiff(1:p, I)] <- Ylabels
                vlabels[I] <- Ilabels
              } else {
                if (is.character(I)) {
                  if (any(is.na(match(I, vlabels))))
                    stop("some vertex labels in 'I' do not form part from the vertex labels in the row and column names from 'g'.")
                  Ilabels <- I
                } else {
                  if (I < 1 || I > p)
                    stop("'I' contains vertices outside the range 1 .. p, with p x p being the dimension of 'g'.")
                  Ilabels <- vlabels[I]
                }
                Ylabels <- setdiff(vlabels, Ilabels)
              }

              from <- vlabels[row(g)[upper.tri(g) & g]]
              to <- vlabels[col(g)[upper.tri(g) & g]]
              df <- rbind(df, data.frame(from=from, to=to, weight=rep(1, length(from)), stringsAsFactors=FALSE))

            } else {
              vlabels <- sort(unique(as.vector(g)))
              if (class(g[1, 1]) == "character") {
                if (!is.character(I))
                  stop("if edges in 'g' are specified by character vertex labels, then vertices in 'I' should be specified as character vertex labels too.")

                Ilabels <- I
                Ylabels <- setdiff(vlabels, Ilabels) ## Ilabels not in vlables correspond to isolated discrete vertices
                df <- rbind(df, data.frame(from=g[, 1], to=g[, 2], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              } else {
                if (is.character(I))
                  stop("if 'I' specifies character vertex labels, then 'g' should contain character vertex labels too.")

                p <- length(unique(c(vlabels, I)))
                Ilabels <- sprintf("I%d", seq(along=I))
                Ylabels <- sprintf("Y%d", seq(along=setdiff(1:p, I)))
                vlabels <- sprintf("X%d", 1:p)
                vlabels[setdiff(1:p, I)] <- Ylabels
                vlabels[I] <- Ilabels
                df <- rbind(df, data.frame(from=vlabels[g[, 1]], to=vlabels[g[, 2]], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              }
            }

            g <- graphBAM(df, nodes=c(Ilabels, Ylabels))
            graph::nodeDataDefaults(g, "type") <- "continuous"
            graph::nodeData(g, Ilabels, "type") <- "discrete"

            rHMgmm(n=n, g=g, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="graphBAM", g="missing"),
          function(n, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g=n, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="missing", g="graphBAM"),
          function(n, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=1L, g, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="numeric", g="graphBAM"),
          function(n=1, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            rHMgmm(n=as.integer(n), g, rho, a, dLevels, tol, verbose)
          })

setMethod("rHMgmm", signature(n="integer", g="graphBAM"),
          function(n=1, g, rho=0.5, a=1, dLevels=2L, tol=0.001, verbose=FALSE) {
            if (is.na(match("type", names(graph::nodeData(g, graph::nodes(g))))))
              stop("if 'g' is a graphBAM object, then its nodes should contain a 'type' attribute specifying whether each of them is 'discrete' or 'continuous'.")

            vtype <- factor(unlist(nodeData(g, nodes(g), "type")))

            pI <- sum(vtype == "discrete")
            pY <- sum(vtype == "continuous")
            Ilabels <- graph::nodes(g)[vtype == "discrete"]
            Ylabels <- graph::nodes(g)[vtype == "continuous"]

            if (pI < 1)
              stop("pI should be equal or larger than 1.")
            if (pY < 1)
              stop("pY should be equal or larger than 1.")

            p <- as.integer(pI + pY)

            labels <- c(Ilabels, Ylabels)

            if (!is.null(labels))
              if (length(labels) != p)
                stop("when set, vertex labels should contain 'pI' + 'pY' variable names.")

            if (rho <= -1/(pY-1) || rho >= 1)
              stop("'rho' should be a real number such that -1/(pY-1) < 'rho' < 1.")

            if (any(dLevels > 2))
              stop("Only binary variables can be used at the moment.")

            if (any(dLevels < 1))
              stop("Discrete variables should have at least two levels specified in argument 'dLevels'.")

            ed <- edges(g)[Ilabels]
            if (any(sapply(ed, function(x, vt) any(vt[x] == "discrete"), vtype)))
              stop("'g' cannot contain edges between vertices of discrete variables.")

            sim <- list()
            for (i in 1:n) {
              sim.sigma <- qpG2Sigma(g=subGraph(Ylabels, g), rho=rho, verbose=verbose)

              sim[[i]] <- HMgmm(g=g, dLevels=dLevels, a=a, rho=rho, sigma=sim.sigma)
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })


## $ accessor operator
setMethod("$", signature(x="HMgmm"),
          function(x, name) {
              switch(name,
                     X=graph::nodes(x@g),
                     I=graph::nodes(x@g)[x@vtype == "discrete"],
                     Y=graph::nodes(x@g)[x@vtype == "continuous"],
                     p=x@pI+x@pY,
                     pI=x@pI,
                     pY=x@pY,
                     g=x@g,
                     mean=function(...) calculateCondMean(x, ...),
                     sigma=x@sigma,
                     a=x@a,
                     eta2=calculateEta2(x),
                     stop("unknown HMgmm slot or parameter. Use names() to find out which are the valid ones.")
                     )
          })

## names method
setMethod("names", signature(x="HMgmm"),
          function(x) {
            c("X", "I", "Y", "p", "pI", "pY", "g", "mean", "sigma", "a", "eta2")
          })

## dim method
setMethod("dim", signature(x="HMgmm"),
          function(x) {
            c(pI=x@pI, pY=x@pY)
          })

## dimnames method
setMethod("dimnames", signature(x="HMgmm"),
          function(x) {
            list(I=x$I, Y=x$Y)
          })

## internal eta2 function
calculateEta2 <- function(object) {
  if (!is.null(object@eta2$value))
    return(object@eta2$value)

    Y <- graph::nodes(object@g)[object@vtype == "continuous"]
    I <- graph::nodes(object@g)[object@vtype == "discrete"]

    eta2 <- rep(NA, object@pY)
    names(eta2) <- Y
    z1xz1 <- YxI <- rep(0, object@pY)
    names(z1xz1) <- names(YxI) <- Y
    z1 <- rep(1, object@pY)
    names(z1) <- Y
    ed <- graph::edges(object@g)
    for (Ii in I) {
      YxIi <- ed[[Ii]] ## assuming discrete variables can only be connected to continuous ones
      nYxIi <- length(YxIi)
      YxI[YxIi] <- YxI[YxIi] + 1

      if (all(YxI <= 1) && nYxIi > 0) {
        z1[YxIi] <- as.vector(solve(object@sigma[YxIi, YxIi, drop=FALSE]) %*% object@a[YxIi])

        if (nYxIi > 1) {
          for (j in 1:(nYxIi-1))
            for (k in (j+1):nYxIi)
              z1xz1 <- z1xz1 + (object@sigma[, YxIi[j]] * object@sigma[, YxIi[k]]) * (z1[YxIi[j]] * z1[YxIi[k]])
        }
      }
    }

    if (all(YxI <= 1)) {
      eta2 <- (object@a^2) / (4*sapply(Y, function(Yi) object@sigma[Yi, Yi] + (1/2)*z1xz1[Yi] + (1/4)*sum(sapply(Y, function(Yj) (z1[Yj] * object@sigma[Yi, Yj])^2))))
      names(eta2) <- Y
    } else
      warning("eta2 can only be computed when every continuous variable is connected in 'g' to at most one discrete variable.")

  object@eta2$value <- eta2
  eta2
}

## internal conditional mean function
calculateCondMean <- function(x, i) {
  if (missing(i)) {
    if (x@pI < 5) {
      i <- 1:(x@dLevels^x$pI)
    } else
      stop("the input HMgmm has too many joint levels to calculate all conditional mean vectors, please specify which are the ones to calculate\n")
  }

  if (class(i) == "numeric" || class(i) == "integer") {
    if (any(i < 1) || any(log(i) > x@pI*log(x@dLevels)))
      stop(sprintf("'i' is either < 1 or larger than %d^%d", x@dLevels, x@pI))

    if (x@dLevels == 2) {
      i <- sapply(i, function(dec) sapply(strsplit(as.character(rev(intToBits(dec-1))),""),`[[`,2)[(32-x@pI+1):32],
                  simplify=FALSE)
      i <- do.call("rbind", i)
    } else {
      i <- sapply(i, function(dec, x) {
                         intf <- dec-1
                         nonbits <- rep("0", x@pI)
                         k <- x@pI
                         while (intf > 0) {
                           nonbits[k] <- as.character(intf %% x@dLevels)
                           intf <- floor(intf / x@dLevels)
                           k <- k - 1
                         }
                         nonbits
                       }, x, simplify=FALSE)
      i <- do.call("rbind", i)
    }
  }

  if (class(i) == "character" || class(i) == "factor") {
    if (length(i) != x@pI)
      stop("'i' should be either a matrix or a data.frame object with discrete joint levels on the rows and discrete variables on the columns.")
    i <- matrix(i, nrow=1, ncol=x@pI, byrow=TRUE)
  }

  if (class(i) != "matrix" && class(i) != "data.frame")
    stop("'i' should be either a matrix or a data.frame object with discrete joint levels on the rows and discrete variables on the columns.")

  p <- x@pI + x@pY
  Y <- graph::nodes(x@g)[x@vtype == "continuous"]
  I <- graph::nodes(x@g)[x@vtype == "discrete"]
  n <- nrow(i)

  if (ncol(i) != x@pI)
    stop("the input discrete data in 'i' should have as many columns as discrete variables in the input HMgmm 'x'.")

  if (is.null(colnames(i)))
    colnames(i) <- I
  else if (!identical(colnames(i), I))
    stop("column names in 'i' should correspond to the discrete variable names in 'x'.")

  i <- matrix(as.integer(factor(i)), ncol=x@pI, dimnames=list(NULL, colnames(i)))

  ## we use the same data matrix mu to store the intermediate values required
  ## during the calculations
  
  mu <- matrix(1, nrow=n, ncol=x@pY, dimnames=list(1:n, Y))

  ## we start by setting on the continuous variables what data points are
  ## associated to a particular discrete level, i.e., from what Gaussian
  ## distribution the mean should be calculated 
  YxI <- Y[which(sapply(graph::edges(x@g, Y), function(xYk, vt) sum(vt[xYk] == "discrete"), x@vtype) > 0)]
  mu[, YxI] <- sapply(YxI, function(Yk, i, pI, mod) {
                               IxYk <- intersect(graph::edges(x@g)[[Yk]], I)
                               apply(i[, IxYk, drop=FALSE], 1,
                                     function(k) sum((k-1)*2^((length(k)-1):0)))
                             }, i, x@pI, x)

  ## calculate canonical parameter h(i) from the homogeneous mixed graphical model
  ## by now this only works with at most one discrete variable associated to a continuous one

  IxY <- I[which(sapply(graph::edges(x@g)[I], function(xIi, vt) sum(vt[xIi] == "continuous"), x@vtype) > 0)]
  for (Ii in IxY) {
    YxIi <- intersect(graph::edges(x@g)[[Ii]], Y)
    if (length(YxIi) > 0) {
      difh <- solve(x@sigma[YxIi, YxIi, drop=FALSE]) %*% x@a[YxIi]
      rownames(difh) <- YxIi

      for (j in YxIi)
        mu[, j] <- mu[, j] * difh[j, ] 
    }
  }

  ## mu = Sigma x h where corresponding h(i) are stored in mu itself for convenience
  mu <- t(as.matrix(x@sigma) %*% t(mu))

  mu
}

## rcmvnorm() uses the rmvnorm() function from the mvtnorm package to sample multivariate normal observations
## with means conditioned in the joint distribution of discrete r.v.'s as defined by an HMgmm
setMethod("rcmvnorm", signature(n="ANY", model="HMgmm"),
          function(n, model, ...) {
            p <- model@pI + model@pY
            Y <- graph::nodes(model@g)[model@vtype == "continuous"]
            I <- graph::nodes(model@g)[model@vtype == "discrete"]
            
            X <- matrix(1, nrow=n, ncol=p, dimnames=list(1:n, c(I, Y)))

            ## simulate discrete data uniformly at random
            X[, I] <- sample(1:model@dLevels, size=n*model@pI, replace=TRUE)

            ## we use the same data matrix X to store the mean values employed
            ## during the simulation process.
            X[, Y] <- calculateCondMean(model, X[, I, drop=FALSE])

            YxI <- Y[which(sapply(graph::edges(model@g)[Y], function(xYi, vt) sum(vt[xYi] == "discrete"), model@vtype) > 0)]
            xtab <- tapply(1:n, apply(X[, YxI, drop=FALSE], 1, function(i) paste(i, collapse="")))
            xtab <- split(as.data.frame(X[, I]), xtab)
            for (i in 1:length(xtab)) {
              li <- xtab[[i]]
              which_n <- as.numeric(rownames(li))
              X[which_n, Y] <- mvtnorm::rmvnorm(length(which_n), mean=X[which_n[1], Y],
                                                sigma=as.matrix(model@sigma), ...)
            }

            X
          })

## show method
setMethod("show", signature(object="HMgmm"),
          function(object) {
            cat(sprintf("\n  Homogeneous mixed graphical Markov model\n  with %d discrete and %d continuous r.v., and %d edges.\n\n",
                        sum(object@vtype == "discrete"), sum(object@vtype == "continuous"), graph::numEdges(object@g)))
            invisible(object)
          })

## summary method
setMethod("summary", signature(object="HMgmm"),
          function(object) {
            ne <- graph::numEdges(object@g)/2
            den <- (100*ne) / choose(object@pI+object@pY, 2)
            nme <- sum(sapply(graph::edges(object@g)[object$I],
                              function(v) length(v) > 0))
            denIxY <- (100*nme) / (object@pI*object@pY)
            adjmY <- as(object@g, "matrix") == 1
            adjmY <- adjmY[object@vtype == "continuous", object@vtype == "continuous"]
            denY <- (100*sum(adjmY)/2) / choose(object@pY, 2)
            deg <- as.integer(graph::degree(object@g))
            macor <- Matrix::cov2cor(object@sigma)[upper.tri(adjmY) & adjmY]
            pacor <- Matrix::cov2cor(Matrix::solve(object@sigma))[upper.tri(adjmY) & adjmY]
            maskYxI <- sapply(graph::nodes(object@g)[object@vtype == "continuous"],
                              function(v, e, vt) any(vt[e[[v]]] == "discrete"),
                              graph::edges(object@g, object$Y), object@vtype)
            a <- object@a[maskYxI]
            new("HMgmmSummary", model=object, density=den, densityIxY=denIxY, densityY=denY,
                degree=deg, macor=macor, pacor=pacor, a=a)
          })


## plot method
setMethod("plot", signature(x="HMgmm"),
          function(x, layoutType="dot", lwd=1, ...) {
            g <- x@g
            g <- Rgraphviz::layoutGraph(g, layoutType=layoutType) ## laying out the graph should happen first otherwise
                                                                  ## it overrrides some of the node parameters below
            graph::nodeRenderInfo(g) <- list(label=do.call("names<-", list(x$X, x$X)),
                                             fill=do.call("names<-",
                                                          list(c(rep("black", sum(x@vtype == "discrete")),
                                                                 rep("white", sum(x@vtype == "continuous"))),
                                                               c(do.call("c", as.list(graph::nodes(g)[x@vtype == "discrete"])),
                                                                 do.call("c", as.list(graph::nodes(g)[x@vtype == "continuous"]))))),
                                             textCol=do.call("names<-",
                                                              list(c(rep("white", sum(x@vtype == "discrete")),
                                                                     rep("black", sum(x@vtype == "continuous"))),
                                                                   c(do.call("c", as.list(graph::nodes(g)[x@vtype == "discrete"])),
                                                                     do.call("c", as.list(graph::nodes(g)[x@vtype == "continuous"]))))))
            mixedEdges <- graph::edges(g, x$I)
            mixedEdges <- paste(rep(names(mixedEdges), times=sapply(mixedEdges, length)),
                                unlist(mixedEdges, use.names=FALSE), sep="~")
            ## it seems Rgraphviz::layoutGraph() alphabetically sorts vertices within edge names so we have to
            ## order the endpoints alphabetically
            mixedEdges <- strsplit(mixedEdges, "~")
            mixedEdges <- lapply(mixedEdges, sort)
            mixedEdges <- sapply(mixedEdges, paste, collapse="~")
            graph::edgeRenderInfo(g) <- list(arrowhead="none", arrowtail="none", lwd=lwd)
            maskhead <- !is.na(match(sapply(strsplit(names(graph::edgeRenderInfo(g)$arrowhead[mixedEdges]), "~"), function(x) x[1]), x$I))
            if (any(maskhead))
              graph::edgeRenderInfo(g) <- list(arrowhead=do.call("names<-",
                                                                 list(rep("open", length(mixedEdges[maskhead])), mixedEdges[maskhead])))
            masktail <- !is.na(match(sapply(strsplit(names(graph::edgeRenderInfo(g)$arrowtail[mixedEdges]), "~"), function(x) x[2]), x$I))
            if (any(masktail))
              graph::edgeRenderInfo(g) <- list(arrowtail=do.call("names<-",
                                                                 list(rep("open", length(mixedEdges[masktail])), mixedEdges[masktail])))
            ## ## it seems Rgraphviz::layoutGraph() alphabetically sorts vertices within edge names so we have to
            ## ## check on what side of the edge names are the discrete variables to decide whether we draw
            ## ## open arrow heads or open arrow tails (sigh!)
            ## if (any(!is.na(match(sapply(strsplit(names(graph::edgeRenderInfo(g)$arrowhead), "~"), function(x) x[1]), x$I)))) {
            ##   graph::edgeRenderInfo(g) <- list(arrowhead="none", arrowtail="none", lwd=lwd)
            ##   graph::edgeRenderInfo(g) <- list(arrowhead=do.call("names<-",
            ##                                                      list(rep("open", length(mixedEdges)), mixedEdges)))
            ## } else {
            ##   graph::edgeRenderInfo(g) <- list(arrowtail="none", arrowtail="none", lwd=lwd)
            ##   graph::edgeRenderInfo(g) <- list(arrowtail=do.call("names<-",
            ##                                                      list(rep("open", length(mixedEdges)), mixedEdges)))
            ## }
            Rgraphviz::renderGraph(g, ...)
         })

## show method for summary
setMethod("show", signature(object="HMgmmSummary"),
          function(object) {
            cat(sprintf("\n  Homogeneous mixed graphical Markov model\n  with %d discrete and %d continuous r.v., and %d edges.\n\n",
                        sum(object@model@vtype == "discrete"), sum(object@model@vtype == "continuous"), graph::numEdges(object@model@g)/2))
            denstr <- ifelse(object@density < 1, sprintf("%.g%%", object@density), sprintf("%.0f%%", object@density))
            denIxYstr <- ifelse(object@densityIxY < 1, sprintf("%.g%%", object@densityIxY), sprintf("%.0f%%", object@densityIxY))
            denYstr <- ifelse(object@densityY < 1, sprintf("%.g%%", object@densityY), sprintf("%.0f%%", object@densityY))
            cat(sprintf("  Graph density: %s (all edges) %s (mixed edges) %s (continuous edges)\n", denstr, denIxYstr, denYstr))
            cat("\n  Degree distribution of the vertices in the graph:\n")
            print(summary(object@degree))
            cat("\n  Distribution of marginal correlations for present continuous edges:\n")
            if (length(object@macor) > 0) print(summary(object@macor)) else cat("NA\n")
            cat("\n  Distribution of partial correlations for present continuous edges:\n")
            if (length(object@pacor) > 0) print(summary(object@pacor)) else cat("NA\n")
            cat("\n  Distribution of additive linear effects for present mixed edges:\n")
            if (length(object@a) > 0) print(summary(object@a)) else cat("NA\n")
            cat("\n")
            invisible(object)
          })
