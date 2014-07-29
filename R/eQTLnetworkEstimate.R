eQTLnetworkEstimate <- function(modelFormula="~ marker", param, firsteQTLnetwork=NULL,
                                verbose=TRUE, BPPARAM=bpparam("SerialParam")) {
  if (!is(param, "eQTLnetworkEstimationParam"))
    stop("argument 'param' must be an 'eQTLnetworkEstimationParam' object.")

  tryCatch({
    rhs <- deparse(formula(modelFormula)[[2]])
  }, error=function(err) {
    cat("Malformed model formula\n")
    stop(conditionMessage(err), call.=FALSE)
  })

  rhs <- strsplit(rhs, " *\\| *")[[1]]
  if (length(rhs) > 2)
    stop("Conditioning operator '|' can only occur once in the model formula.")

  rhs <- strsplit(rhs, " *\\+ *")
  if (length(rhs) == 1)
    rhs[[2]] <- character(0)

  names(rhs) <- c("explanatory" ,"conditioning")
  ## we'll need this when we handle interactions in the future
  ## rhs[[1]] <- unlist(lapply(rhs[[1]], strsplit, " *\\* *| *: *"), recursive=FALSE)
  ## rhs[[2]] <- unlist(lapply(rhs[[2]], strsplit, " *\\* *| *: *"), recursive=FALSE)

  clusterSize <- 1 ## this should change to use BiocParallel
  if (!is(BPPARAM, "SerialParam")) {
    if (!is(BPPARAM, "SnowParam"))
      stop("At the moment parallelization only works with 'SnowParam' parallel backends.")

    if (BPPARAM$.clusterargs$type != "MPI")
      stop("At the moment the only type of parallel cluster available is 'MPI'") 

    if (BPPARAM$.clusterargs$spec > 0)
      clusterSize <- BPPARAM$.clusterargs$spec
  }

  explanatory <- unlist(lapply(.expandTerms(rhs$explanatory, param), function(x) x$v))

  pvaluesG0 <- nrr <- as(Matrix(numeric(), nrow=0, ncol=0), "dspMatrix")

  if (is(firsteQTLnetwork, "eQTLnetwork")) {
    if (!identical(explanatory, unlist(lapply(.expandTerms(firsteQTLnetwork@rhs$explanatory, param),
                                              function(x) x$v))))
      stop("Explanatory terms in the model formula do not match those of the given eQTLnetwork estimate\n")
    pvaluesG0 <- firsteQTLnetwork@pvaluesG0
    nrr <- firsteQTLnetwork@nrr
  }

  res <- NULL
  if (length(rhs$conditioning) == 0) {
    if (is(firsteQTLnetwork, "eQTLnetwork")) {
      if (nrow(pvaluesG0) > 0)
        stop("The provided eQTLnetwork has already a G^0 estimate\n")
      rhs$conditioning <- firsteQTLnetwork@rhs@conditioning
    }

    message("Calculating pairwise (conditional) independence tests\n")
    all.marginal.tests <- qpAllCItests(ggData(param), I=param@dVars, Q=NULL,
                                       pairup.i=geneNames(param), pairup.j=explanatory,
                                       return.type="p.value", clusterSize=clusterSize,
                                       verbose=verbose)

    res <- new("eQTLnetwork", geneticMap=geneticMap(param), physicalMap=physicalMap(param),
               organism=param@organism, genome=param@genome, geneAnnotation=geneAnnotation(param),
               geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
               pvaluesG0=all.marginal.tests$p.value, nrr=nrr, modelFormula=modelFormula,
               rhs=rhs, qOrders=integer())
  } else {
    eterms <- .expandTerms(rhs$conditioning, param)
    maskq0 <- sapply(eterms, function(x) all(x$q==0))
    fix.Q <- unlist(lapply(eterms[maskq0], '[[', "v"))
    restrict.Q <- unlist(lapply(eterms[!maskq0], '[[', "v"))
    qorders <- lapply(eterms[!maskq0], function(x) x$q)

    ## check that q-orders associated to different terms are identical
    if (!all(sapply(qorders, function(x, r) identical(x, r), qorders[[1]])))
      stop("q values should be identical throughout the given conditioning terms")
    qorders <- sort(qorders[[1]])

    if (is(firsteQTLnetwork, "eQTLnetwork")) {
      if (length(intersect(qorders, firsteQTLnetwork@qOrders)) > 0)
        stop("The provided eQTLnetwork was estimated using already some of the given q-orders\n")
    }

    if (any(qorders == 0)) {
      if (is(firsteQTLnetwork, "eQTLnetwork")) {
        if (nrow(pvaluesG0) > 0)
          stop("The provided eQTLnetwork has already a G^0 estimate\n")
      }
      message("Calculating pairwise (conditional) independence tests\n")
      pvaluesG0 <- qpAllCItests(ggData(param), I=param@dVars, Q=NULL,
                                pairup.i=geneNames(param), pairup.j=explanatory,
                                return.type="p.value", clusterSize=clusterSize,
                                verbose=verbose)$p.value
      qorders <- qorders[qorders > 0]
    }
    k <- 0
    if (is(firsteQTLnetwork, "eQTLnetwork"))
      k <- length(firsteQTLnetwork@qOrders)

    for (i in seq(along=qorders)) {
      q <- qorders[i]
      message(sprintf("Estimating non-rejection rates with q=%d\n", q))
      nrr.q <- qpNrr(ggData(param), I=param@dVars, q=q, pairup.i=geneNames(param), pairup.j=explanatory,
                     restrict.Q=restrict.Q, fix.Q=fix.Q, clusterSize=clusterSize, verbose=verbose)
      if (i+k == 1)
        nrr <- nrr.q
      else
        nrr <- (nrr.q + nrr*(i+k-1)) / (i+k)
    }
    if (is(firsteQTLnetwork, "eQTLnetwork"))
      qorders <- sort(c(qorders, firsteQTLnetwork@qOrders))

    res <- new("eQTLnetwork", geneticMap=geneticMap(param), physicalMap=physicalMap(param),
               organism=param@organism, genome=param@genome, geneAnnotation=geneAnnotation(param),
               geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
               pvaluesG0=pvaluesG0, nrr=nrr, modelFormula=modelFormula, rhs=rhs, qOrders=qorders)
  }

  res
}

.expandTerms <- function(termvec, param) {
  varsnqs <- lapply(termvec,
                    function(x, param) {
                      y <- list(v=x, q=0)
                      m <- regexec("([^\\(]+)\\(q = ([0-9, ]+)\\)", x)
                      m <- regmatches(x, m)[[1]]
                      if (length(m) > 0) {
                        y$v <- m[2]
                        y$q <- as.integer(strsplit(m[3], ",")[[1]])
                      }

                      if (y$v == "marker")
                        y$v <- markerNames(param)
                      else if (y$v == "gene")
                        y$v <- geneNames(param)
                      else if (y$v %in% colnames(mcols(geneAnnotation(param)))) {
                        if (!is(geneAnnotation(param)[[y$v]], "logical"))
                            stop(sprintf("%s is not a column of logical values in the gene annotation", y$v))
                        y$v <- names(geneAnnotation(param))[geneAnnotation(param)[[y$v]]]
                      }
                      y
                    }, param)
  varsnqs
}
