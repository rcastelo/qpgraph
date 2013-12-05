## graphParam constructors

erGraphParam <- function(p=4L, m=4L, prob=NA_real_, labels=as.character(1:p)) {
  new("erGraphParam", p=as.integer(p), m=as.integer(m), prob=prob, labels=labels)
}

dRegularGraphParam <- function(p=4L, d=2L, exclude=as.integer(NULL), labels=as.character(1:p)) {
  if ((p*d) %% 2 != 0)
    stop("The number of vertices p times the degree d of each vertex, i.e., the product p x d, should be even in order to sample a d-regular graph on p vertices uniformly at random\n")

  new("dRegularGraphParam", p=as.integer(p), d=as.integer(d), exclude=exclude, labels=labels)
}

erMarkedGraphParam <- function(pI=1L, pY=3L, m=4L, prob=NA_real_,
                               Ilabels=paste0("I", 1:pI), Ylabels=paste0("Y", 1:pY)) {
  new("erMarkedGraphParam", p=as.integer(pI+pY), pI=as.integer(pI), pY=as.integer(pY), m=as.integer(m),
      prob=prob, Ilabels=Ilabels, Ylabels=Ylabels, labels=c(Ilabels, Ylabels))
}

dRegularMarkedGraphParam <- function(pI=1L, pY=3L, d=2L, exclude=as.integer(NULL),
                                     Ilabels=paste0("I", 1:pI), Ylabels=paste0("Y", 1:pY)) {
  if (((pI+pY)*d) %% 2 != 0)
    stop("The number of vertices p times the degree d of each vertex, i.e., the product p x d, should be even in order to sample a d-regular graph on p vertices uniformly at random\n")

  new("dRegularMarkedGraphParam", p=as.integer(pI+pY), pI=as.integer(pI), pY=as.integer(pY),
      d=as.integer(d), exclude=exclude, Ilabels=Ilabels, Ylabels=Ylabels, labels=c(Ilabels, Ylabels))
}


## graphParam show method
setMethod("show", signature(object="graphParam"),
          function(object) {
            graphtype <- "pure"
            if (is(object, "markedGraphParam"))
              graphtype <- "marked"
            graphmodel <- "Erdos-Renyi"
            if (is(object, "dRegularGraphParam") || is(object, "dRegularMarkedGraphParam"))
              graphmodel <- "d-regular"

            cat(sprintf("\n  %s %s graph parameter object\n", graphmodel, graphtype))
            cat(sprintf("  No. of %s vertices: %d\n", graphtype, object@p))

            if (graphtype == "marked") {
              cat(sprintf("  No. of dot (I) vertices: %d\n", object@pI))
              cat(sprintf("  No. of circle (Y) vertices: %d\n", object@pY))
            }

            if (graphmodel == "Erdos-Renyi") {
              if (is.na(object@prob))
                cat(sprintf("  No. of edges: %d\n", object@m))
              else
                cat(sprintf("  Edge probability: %.3f\n", object@prob))
            } else if (graphmodel == "d-regular") {
              cat(sprintf("  Constant degree: %d\n", object@d))
              if (length(object@exclude) > 0)
                cat(sprintf("  No. of excluded vertices: %d\n", length(object@exclude)))
            } else stop("unknown graph model.")

            if (graphtype == "pure") {
              vtcstr <- paste(object@labels, collapse=", ")
              if (object@p > 6)
                vtcstr <- sprintf("%s ...", paste(head(object@labels), collapse=", "))
              cat(sprintf("  Vertex labels: %s\n\n", vtcstr))
            } else {
              vtcstr <- paste(object@Ilabels, collapse=", ")
              if (object@pI > 6)
                vtcstr <- sprintf("%s ...", paste(head(object@Ilabels), collapse=", "))
              cat(sprintf("  Dot (I) vertex labels: %s\n", vtcstr))
              vtcstr <- paste(object@Ylabels, collapse=", ")
              if (object@pY > 6)
                vtcstr <- sprintf("%s ...", paste(head(object@Ylabels), collapse=", "))
              cat(sprintf("  Circle (Y) vertex labels: %s\n\n", vtcstr))
            }

            invisible(object)
          })

## graphParam plot method
setMethod("plot", signature(x="graphBAM"),
          function(x, layoutType="dot", ...) {
            if (!is.na(match("type", names(nodeData(x))))) {
              vtype <- unlist(graph::nodeData(x, graph::nodes(x), "type"))
              vlabels <- graph::nodes(x)

              g <- x
              g <- Rgraphviz::layoutGraph(g, layoutType=layoutType)
              graph::nodeRenderInfo(g) <- list(label=do.call("names<-", list(vlabels, vlabels)),
                                               fill=do.call("names<-",
                                                            list(c(rep("black", sum(vtype == "discrete")),
                                                                   rep("white", sum(vtype == "continuous"))),
                                                                 c(do.call("c", as.list(vlabels[vtype == "discrete"])),
                                                                   do.call("c", as.list(vlabels[vtype == "continuous"]))))),
                                               textCol=do.call("names<-",
                                                            list(c(rep("white", sum(vtype == "discrete")),
                                                                   rep("black", sum(vtype == "continuous"))),
                                                                 c(do.call("c", as.list(vlabels[vtype == "discrete"])),
                                                                   do.call("c", as.list(vlabels[vtype == "continuous"]))))))
              Rgraphviz::renderGraph(g, ...)
            } else
              invisible(callNextMethod())
          })

## graph simulation methods
setMethod("rgraphBAM", signature(n="graphParam", param="missing"),
          function(n, param) {
            rgraphBAM(n=1L, n)
          })

setMethod("rgraphBAM", signature(n="missing", param="graphParam"),
          function(n, param) {
            rgraphBAM(n=1L, param)
          })

setMethod("rgraphBAM", signature(n="numeric", param="graphParam"),
          function(n=1, param) {
            rgraphBAM(as.integer(n), param)
          })

setMethod("rgraphBAM", signature(n="integer", param="erGraphParam"),
          function(n=1L, param) {
            p <- param@p
            m <- param@m
            prob <- param@prob
            labels <- param@labels

            sim <- list()
            for (i in 1:n) {
              g <- NA
              if (is.na(prob))
                g <- graph::randomEGraph(labels, edges=m)
              else
                g <- graph::randomEGraph(labels, p=prob)

              g <- as(g, "graphBAM")

              if (is(param, "markedGraphParam")) {
                nodeDataDefaults(g, "type") <- "continuous"
                nodeData(g, param@Ilabels, "type") <- "discrete"
              }

              sim[[i]] <- as(g, "graphBAM")
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

setMethod("rgraphBAM", signature(n="integer", param="dRegularGraphParam"),
          function(n=1L, param, verbose=FALSE, R.code.only=FALSE) {
            p <- param@p
            d <- param@d
            exclude <- param@exclude
            labels <- param@labels

            sim <- list()
            for (i in 1:n) {
              g <- qpgraph:::qpRndRegularGraph(p, d, labels, exclude, verbose, return.type="graphBAM", R.code.only)
              if (is(param, "markedGraphParam")) {
                nodeDataDefaults(g, "type") <- "continuous"
                nodeData(g, param@Ilabels, "type") <- "discrete"
              }

              sim[[i]] <- g
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })
