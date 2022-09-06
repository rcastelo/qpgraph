## constructor methods
setMethod("qpGraph", signature(nrrMatrix="matrix"),
          function(nrrMatrix, ...) {
            otherArgs <- names(list(...))
            if ("threshold" %in% otherArgs || "return.type" %in% otherArgs)
              return(.old_qpGraph(nrrMatrix, ...))

            .qpGraph(as(as(as(nrrMatrix, "dMatrix"), "symmetrixMatrix"), "packedMatrix"), ...)
          })

setMethod("qpGraph", signature(nrrMatrix="dspMatrix"),
          function(nrrMatrix, ...) {
            otherArgs <- names(list(...))
            if ("threshold" %in% otherArgs || "return.type" %in% otherArgs)
              return(.old_qpGraph(nrrMatrix, ...))

            .qpGraph(nrrMatrix, ...)
          })

## names method
setMethod("names", signature(x="qpGraph"),
          function(x) {
            c("X", "p", "q", "n", "epsilon", "g")
          })

## $ accessor operator
setMethod("$", signature(x="qpGraph"),
          function(x, name) {
            switch(name,
                   X=graph::nodes(x@g),
                   p=x@p,
                   q=x@q,
                   n=x@n,
                   epsilon=x@epsilon,
                   g=x@g,
                   stop("unknown qpGraph slot or parameter. Use names() to find out which are the valid ones.")
                   )
          })

## mapping of non-negative integers to unordered edges of non-negative integer vertices
.i2e <- function(i) {
  v <- 1 + floor(-0.5+sqrt(0.25+2*i))
  w <- i - (v*(v-1))/2
  cbind(v, w)
}

.qpGraph <- function(nrrMatrix, epsilon=NA_real_, topPairs=NA_integer_,
                     pairup.i=NULL, pairup.j=NULL, q, n) {
  p <- nrow(nrrMatrix)

  if (!isSymmetric(nrrMatrix))
    stop("'nrrMatrix' is not symmetric.")

  if (missing(n))
    n <- NA_integer_
  if (missing(q))
    q <- NA_integer_

  vertex.labels <- NULL
  if (is.null(colnames(nrrMatrix))) {
    vertex.labels <- as.character(1:p)
  } else {
    vertex.labels <- colnames(nrrMatrix)
  }

  if (is.na(epsilon) && is.na(topPairs))
    stop("either 'epsilon' or 'topPairs' should be set different to NULL\n")

  if (!is.na(epsilon) && !is.na(topPairs))
    stop("only either 'epsilon' or 'topPairs' can be set different to NULL\n")

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
       (is.null(pairup.i) && !is.null(pairup.j)))
    stop("'pairup.i' and 'pairup.j' should both either be set to NULL or contain subsets of variables\n")

  if (!is.null(pairup.i) && !is.null(pairup.j))  {
    if (is.null(colnames(nrrMatrix)))
      stop("when using 'pairup.i' and 'pairup.j', nrrMatrix should have row and column names\n")

    var.names <- colnames(nrrMatrix)
    pairup.i <- match(pairup.i, var.names)
    if (sum(is.na(pairup.i)) > 0)
      stop("'pairup.i' is not a subset of the variables forming the data\n")
    pairup.j <- match(pairup.j, var.names)
    if (sum(is.na(pairup.j)) > 0)
      stop("'pairup.j' is not a subset of the variables forming the data\n")

    pairup.ij.int <- intersect(pairup.i, pairup.j)
    pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
    pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

    nomeasurementsMask <- matrix(FALSE,nrow=p,ncol=p)
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

  nrrUT <- nrrMatrix[upper.tri(nrrMatrix)]
  df <- NULL
  if (!is.na(epsilon)) {                            ## epsilon (cutoff on the non-rejection rate)
    idx <- which(nrrUT <= epsilon)
    idx <- .i2e(idx-1) + 1
    df <- data.frame(from=vertex.labels[idx[, 1]],
                     to=vertex.labels[idx[, 2]],
                     weight=rep(1, nrow(idx)))
  } else {                                          ## topPairs
    nrrUTsorted <- sort(nrrUT, partial=topPairs)[1:topPairs]
    idx <- which(nrrUT %in% nrrUTsorted)
    if (length(idx) > topPairs)
      idx <- idx[order(nrrUT[idx])][1:topPairs]     ## handle when two NRR values are identical
    idx <- .i2e(idx-1) + 1
    df <- data.frame(from=vertex.labels[idx[, 1]],
                     to=vertex.labels[idx[, 2]],
                     weight=rep(1, topPairs))
  }

  g <- graphBAM(df, nodes=vertex.labels)

  qpg <- new("qpGraph", p=as.integer(p), n=as.integer(n),
             q=as.integer(q), epsilon=as.numeric(epsilon),
             g=g)
  qpg
}

## show method
setMethod("show", signature(object="qpGraph"),
          function(object) {
            cat(sprintf("\n qp-graph object with p=%d and q={%s}\n\n", object@p, paste(object@q, collapse=", ")))
            invisible(object)
          })
