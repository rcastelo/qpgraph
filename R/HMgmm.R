setValidity("HMgmm",
            function(object) {
              valid <- TRUE

              if (graph::edgemode(object@g) != "undirected")
                valid <- "currently only undirected homgenous mixed graphical Markov models are supported and 'g' is not undirected."

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

            a <- rep(0, pY)
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
            if (length(a) == 1)
              a <- rep(a, pY)
            if (is.null(names(a)))
              names(a) <- nodes(g)[vtype == "continuous"]
            new("HMgmm", pI=pI, pY=pY, g=g, vtype=vtype, dLevels=dLevels, a=a,
                rho=rho, sigma=as(sigma, "dspMatrix"), eta2=new.env(parent=emptyenv()))
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
  YxI <- Y[which(sapply(graph::edges(x@g)[Y], function(xYk, vt) sum(vt[xYk] == "discrete"), x@vtype) > 0)]
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
            ne <- graph::numEdges(object@g)
            den <- (100*ne) / choose(object@pI+object@pY, 2)
            nme <- sum(sapply(graph::edges(object@g)[object@vtype == "discrete"],
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
                              graph::edges(object@g), object@vtype)
            a <- object@a[maskYxI]
            new("HMgmmSummary", model=object, density=den, densityIxY=denIxY, densityY=denY,
                degree=deg, macor=macor, pacor=pacor, a=a)
          })


## plot method
setMethod("plot", signature(x="HMgmm"),
          function(x, layoutType="dot") {
            g <- x@g
            g <- Rgraphviz::layoutGraph(g, layoutType=layoutType) ## this should come first otherwise it overrrides some of the node parameters below
            graph::nodeRenderInfo(g) <- list(shape="ellipse", fixedsize=FALSE,
                                             label=do.call("names<-", list(x$X, x$X)),
                                             fill=do.call("names<-",
                                                              list(c(rep("black", sum(x@vtype == "discrete")), rep("white", sum(x@vtype == "continuous"))),
                                                                   c(do.call("c", as.list(nodes(g)[x@vtype == "discrete"])), do.call("c", as.list(nodes(g)[x@vtype == "continuous"]))))),
                                             textCol=do.call("names<-",
                                                              list(c(rep("white", sum(x@vtype == "discrete")), rep("black", sum(x@vtype == "continuous"))),
                                                                   c(do.call("c", as.list(nodes(g)[x@vtype == "discrete"])), do.call("c", as.list(nodes(g)[x@vtype == "continuous"]))))))
            Rgraphviz::renderGraph(g)
         })

## show method for summary
setMethod("show", signature(object="HMgmmSummary"),
          function(object) {
            cat(sprintf("\n  Homogeneous mixed graphical Markov model\n  with %d discrete and %d continuous r.v., and %d edges.\n\n",
                        sum(object@model@vtype == "discrete"), sum(object@model@vtype == "continuous"), graph::numEdges(object@model@g)))
            denstr <- ifelse(object@density < 1, sprintf("%.g%%", object@density), sprintf("%.0f%%", object@density))
            denIxYstr <- ifelse(object@densityIxY < 1, sprintf("%.g%%", object@densityIxY), sprintf("%.0f%%", object@densityIxY))
            denYstr <- ifelse(object@densityY < 1, sprintf("%.g%%", object@densityY), sprintf("%.0f%%", object@densityY))
            cat(sprintf("  Graph density: %s (all edges) %s (mixed edges) %s (continuous edges)\n", denstr, denIxYstr, denYstr))
            cat("\n  Degree distribution of the undirected graph:\n")
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
