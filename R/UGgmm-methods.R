setValidity("UGgmm",
            function(object) {
              valid <- TRUE

              if (graph::edgemode(object@g) != "undirected")
                valid <- "currently only undirected Gaussian graphical models are supported and 'g' is not undirected."

              if (class(valid) == "logical" && length(object@mean) != object@p)
                valid <- "'mean' should contain 'p' real numbers corresponding to the mean vector parameter."

              if (class(valid) == "logical" && graph::numNodes(object@g) != object@p)
                valid <- "the number of vertices in 'g' should equal the length of the 'mean' vector."

              p <- (d <- dim(object@sigma))[1]
              if (class(valid) == "logical" && p != d[2])
                valid <- "'sigma' should be a squared matrix."

              if (class(valid) == "logical" && p != object@p)
                valid <- "the dimension of sigma should match the length of the mean vector."

              if (class(valid) == "logical" && !isSymmetric(as(object@sigma, "matrix")))
                valid <- "'sigma' should be symmetric."

              valid
            })
 
## constructor methods
setMethod("UGgmm", signature(g="missing"),
          function(g, verbose) {
            g <- graphBAM(as.data.frame(matrix(NA, nrow=0, ncol=3,
                                               dimnames=list(NULL, c("from", "to", "weight")))),
                          nodes=sprintf("%d", 1:5))
            mean <- rep(0, graph::numNodes(g))
            sigma <- diag(graph::numNodes(g))
            
            UGgmm(g, mean, sigma, verbose)
          })

setMethod("UGgmm", signature(g="matrix"),
          function(g, mean, sigma, verbose) {
            p <- (d <- dim(g))[1]
            if (p != d[2] && d[2] != 2)
              stop("If 'g' is a matrix it should be either a squared symmetric adjacency matrix or a two-column matrix with un-ordered vertex pairs in the rows defining the edge set of an undirected graph.")

            df <- as.data.frame(matrix(NA, nrow=0, ncol=3, dimnames=list(NULL, c("from", "to", "weight"))),
                                stringsAsFactors=FALSE)
            vlabels <- NULL
            if (p == d[2]) {
              if (!isSymmetric(g))
                stop("'g' is not a symmetric matrix\n")

              vlabels <- colnames(g)
              if (is.null(vlabels))
                vlabels <- sprintf("%d", 1:p)
              from <- vlabels[row(g)[upper.tri(g) & g]]
              to <- vlabels[col(g)[upper.tri(g) & g]]
              df <- rbind(df, data.frame(from=from, to=to, weight=rep(1, length(from)), stringsAsFactors=FALSE))
            } else {
              if (class(g[1, 1]) == "character")
                df <- rbind(df, data.frame(from=g[, 1], to=g[, 2], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              else {
                vlabels <- sprintf("%d", sort(unique(as.vector(g))))
                df <- rbind(df, data.frame(from=vlabels[g[, 1]], to=vlabels[g[, 2]], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              }
            }

            UGgmm(graphBAM(df), mean, sigma, verbose)
          })

setMethod("UGgmm", signature(g="graphBAM"),
          function(g, mean=rep(0, graph::numNodes(g)),
                   sigma=diag(graph::numNodes(g)), verbose=FALSE) {
            new("UGgmm", p=length(mean), g=g, mean=mean,
                sigma=as(as(as(sigma, "dMatrix"), "symmetricMatrix"), "packedMatrix"))
          })

## constructor simulation methods
setMethod("rUGgmm", signature(n="graphParam", g="missing"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g=n, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="missing", g="graphParam"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="numeric", g="graphParam"),
          function(n=1, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(as.integer(n), g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="integer", g="graphParam"),
          function(n=1L, g, rho=0.5, tol=0.001, verbose=FALSE) {
            p <- g@p

            if (rho <= -1/(p-1) || rho >= 1)
              stop("'rho' should be a real number such that -1/(p-1) < 'rho' < 1.")

            sim <- list()
            for (i in 1:n) {
              sim.g <- rgraphBAM(g)

              sim.sigma <- qpG2Sigma(sim.g, rho=rho, tol=tol, verbose=verbose)

              sim[[i]] <- UGgmm(g=sim.g, mean=rep(0, p), sigma=sim.sigma)
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

setMethod("rUGgmm", signature(n="matrix", g="missing"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g=n, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="missing", g="matrix"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="numeric", g="matrix"),
          function(n=1, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(as.integer(n), g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="integer", g="matrix"),
          function(n=1L, g, rho=0.5, tol=0.001, verbose=FALSE) {
            p <- (d <- dim(g))[1]
            if (p != d[2] && d[2] != 2)
              stop("If 'g' is a matrix it should be either a squared symmetric adjacency matrix or a two-column matrix with un-ordered vertex pairs in the rows defining the edge set of an undirected graph.")

            df <- as.data.frame(matrix(NA, nrow=0, ncol=3, dimnames=list(NULL, c("from", "to", "weight"))),
                                stringsAsFactors=FALSE)
            vlabels <- NULL
            if (p == d[2]) {
              if (!isSymmetric(g))
                stop("'g' is not a symmetric matrix\n")

              if (class(g[1, 1]) == "integer" || class(g[1, 1]) == "numeric") {
                if (verbose)
                  warning("coercing input numeric adjacency matrix 'g' to a logical adjacency matrix\n")

                g <- g != 0
              }

              vlabels <- colnames(g)
              if (is.null(vlabels))
                vlabels <- sprintf("%d", 1:p)
              from <- vlabels[row(g)[upper.tri(g) & g]]
              to <- vlabels[col(g)[upper.tri(g) & g]]
              df <- rbind(df, data.frame(from=from, to=to, weight=rep(1, length(from)), stringsAsFactors=FALSE))
            } else {
              if (class(g[1, 1]) == "character")
                df <- rbind(df, data.frame(from=g[, 1], to=g[, 2], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              else {
                vlabels <- sprintf("%d", sort(unique(as.vector(g))))
                df <- rbind(df, data.frame(from=vlabels[g[, 1]], to=vlabels[g[, 2]], weight=rep(1, nrow(g)), stringsAsFactors=FALSE))
              }
            }

            rUGgmm(n=n, g=graphBAM(df, nodes=vlabels), rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="graphBAM", g="missing"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g=n, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="missing", g="graphBAM"),
          function(n, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(n=1L, g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="numeric", g="graphBAM"),
          function(n=1, g, rho=0.5, tol=0.001, verbose=FALSE) {
            rUGgmm(as.integer(n), g, rho, tol, verbose)
          })

setMethod("rUGgmm", signature(n="integer", g="graphBAM"),
          function(n=1L, g, rho=0.5, tol=0.001, verbose=FALSE) {
            p <- graph::numNodes(g)

            if (rho <= -1/(p-1) || rho >= 1)
              stop("'rho' should be a real number such that -1/(p-1) < 'rho' < 1.")

            sim <- list()
            for (i in 1:n) {
              sim.sigma <- qpG2Sigma(g, rho=rho, tol=tol, verbose=verbose)

              sim[[i]] <- UGgmm(g=g, mean=rep(0, p), sigma=sim.sigma)
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

## names method
setMethod("names", signature(x="UGgmm"),
          function(x) {
            c("X", "p", "g", "mean", "sigma")
          })

## $ accessor operator
setMethod("$", signature(x="UGgmm"),
          function(x, name) {
            switch(name,
                   X=graph::nodes(x@g),
                   p=x@p,
                   g=x@g,
                   mean=x@mean,
                   sigma=x@sigma,
                   stop("unknown UGgmm slot or parameter. Use names() to find out which are the valid ones.")
                   )
          })

## dim method
setMethod("dim", signature(x="UGgmm"),
          function(x) {
            x@p 
          })

## dimnames method
setMethod("dimnames", signature(x="UGgmm"),
          function(x) {
            graph::nodes(x@g) 
          })

## show method
setMethod("show", signature(object="UGgmm"),
          function(object) {
            cat(sprintf("\n  Undirected Gaussian graphical Markov model\n  with %d r.v. and %d edges.\n\n",
                        graph::numNodes(object@g), graph::numEdges(object@g)))
            invisible(object)
          })

## summary method
setMethod("summary", signature(object="UGgmm"),
          function(object) {
            ne <- graph::numEdges(object@g)
            den <- (100*ne) / choose(object@p, 2)
            deg <- as.integer(graph::degree(object@g))
            ed <- matrix(match(unlist(graph::edges(object@g), use.names=FALSE), graph::nodes(object@g)), ncol=2, byrow=TRUE)
            adjm <- as(object@g, "matrix") == 1
            macor <- Matrix::cov2cor(object@sigma)[upper.tri(adjm) & adjm]
            pacor <- Matrix::cov2cor(Matrix::solve(object@sigma))[upper.tri(adjm) & adjm]
            new("UGgmmSummary", model=object, density=den, degree=deg, macor=macor, pacor=pacor)
          })


## plot method
setMethod("plot", signature(x="UGgmm"),
          function(x, ...) {
            Rgraphviz::plot(x@g, ...)
         })

## overload rmvnorm() from the mvtnorm package to sample multivariate normal observations from an UGgmm
rmvnorm <- function(n, mean, ...) UseMethod("rmvnorm", mean)
rmvnorm.default <- function(n, mean, ...) mvtnorm::rmvnorm(n, mean, ...)
rmvnorm.numeric <- function(n, mean, ...) mvtnorm::rmvnorm(n, mean, ...)
setMethod("rmvnorm", c(n="numeric", mean="numeric"), rmvnorm.numeric)
setMethod("rmvnorm", c(n="integer", mean="numeric"), rmvnorm.numeric)

rmvnorm.UGgmm <- function(n, mean, ...) {
                   X <- mvtnorm::rmvnorm(n=n, mean=mean$mean, sigma=as.matrix(mean$sigma), ...)
                   colnames(X) <- graph::nodes(mean$g)
                   X
                 }
setMethod("rmvnorm", signature=c(n="numeric", mean="UGgmm"), rmvnorm.UGgmm)

## show method for summary
setMethod("show", signature(object="UGgmmSummary"),
          function(object) {
            cat(sprintf("\n  Undirected Gaussian graphical Markov model\n  with %d r.v. and %d edges.\n\n",
                        object@model@p, graph::numEdges(object@model@g)))
            denstr <- ifelse(object@density < 1, sprintf("%.g%%", object@density), sprintf("%.0f%%", object@density))
            cat(sprintf("  Graph density: %s\n", denstr))
            cat("\n  Degree distribution of the undirected graph:\n")
            print(summary(object@degree))
            cat("\n  Distribution of marginal correlations for present edges:\n")
            if (length(object@macor) > 0) print(summary(object@macor)) else cat("NA\n")
            cat("\n  Distribution of partial correlations for present edges:\n")
            if (length(object@pacor) > 0) print(summary(object@pacor)) else cat("NA\n")
            cat("\n")
            invisible(object)
          })
