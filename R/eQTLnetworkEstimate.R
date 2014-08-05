## methods for estimating eQTL networks

setMethod("eQTLnetworkEstimate", signature=c(param="eQTLnetworkEstimationParam",
                                             model="formula",
                                             estimate="missing"),
          function(param, model, estimate,
                   verbose=TRUE, BPPARAM=bpparam("SerialParam")) {
            if (!is(param, "eQTLnetworkEstimationParam"))
              stop("argument 'param' must be an 'eQTLnetworkEstimationParam' object.")

            clusterSize <- 1 ## this should change to use BiocParallel
            if (!is(BPPARAM, "SerialParam")) {
              if (!is(BPPARAM, "SnowParam"))
                stop("At the moment parallelization only works with 'SnowParam' parallel backends.")

              if (BPPARAM$.clusterargs$type != "MPI")
                stop("At the moment the only type of parallel cluster available is 'MPI'") 

              if (BPPARAM$.clusterargs$spec > 0)
                clusterSize <- BPPARAM$.clusterargs$spec
            }

            rhs <- .parseFormula(model)

            pvaluesG0 <- nrr <- as(Matrix(numeric(), nrow=0, ncol=0), "dspMatrix")
            fix.Q <- restrict.Q <- NULL
            qorders <- integer()
            if (length(rhs$conditioning) > 0) {
              maskq0 <- sapply(rhs$conditioning, function(x) all(x$q==0))
              fix.Q <- unlist(lapply(rhs$conditioning[maskq0], '[[', "v"))
              restrict.Q <- unlist(lapply(rhs$conditioning[!maskq0], '[[', "v"))
              qorders <- lapply(rhs$conditioning[!maskq0], function(x) x$q)

              ## check that q-orders associated to different terms are identical
              if (!all(sapply(qorders, function(x, r) identical(x, r), qorders[[1]])))
                stop("q values should be identical throughout the given conditioning terms")
              qorders <- sort(qorders[1][[1]]) ## the [1][[1]] notation handles list() and [[1]] alone does not

              if (is.null(qorders))
                qorders <- 0L
            }

            res <- NULL

            if (length(rhs$conditioning) == 0 || any(qorders == 0)) {
              message("Calculating pairwise (conditional) independence tests\n")
              pvaluesG0 <- qpAllCItests(ggData(param), I=param@dVars, Q=fix.Q,
                                        pairup.i=geneNames(param), pairup.j=rhs$explanatory,
                                        return.type="p.value", clusterSize=clusterSize,
                                        verbose=verbose)$p.value
              qorders <- qorders[qorders > 0]
            }
            
            if (length(rhs$conditioning) > 0) {
              k <- 0

              for (i in seq(along=qorders)) {
                q <- qorders[i]
                message(sprintf("Estimating non-rejection rates with q=%d\n", q))
                nrr.q <- qpNrr(ggData(param), I=param@dVars, q=q, pairup.i=geneNames(param), pairup.j=rhs$explanatory,
                               restrict.Q=restrict.Q, fix.Q=fix.Q, clusterSize=clusterSize, verbose=verbose)
                if (i+k == 1)
                  nrr <- nrr.q
                else
                  nrr <- (nrr.q + nrr*(i+k-1)) / (i+k)
              }
            }

            new("eQTLnetwork", geneticMap=geneticMap(param), physicalMap=physicalMap(param),
                organism=param@organism, genome=param@genome, geneAnnotation=geneAnnotation(param),
                geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
                pvaluesG0=pvaluesG0, nrr=nrr, modelFormula=model, rhs=rhs, qOrders=qorders)
          })

setMethod("eQTLnetworkEstimate", signature=c(param="eQTLnetworkEstimationParam",
                                             model="formula",
                                             estimate="eQTLnetwork"),
          function(param, model, estimate,
                   verbose=TRUE, BPPARAM=bpparam("SerialParam")) {
            if (!is(param, "eQTLnetworkEstimationParam"))
              stop("argument 'param' must be an 'eQTLnetworkEstimationParam' object.")

            clusterSize <- 1 ## this should change to use BiocParallel
            if (!is(BPPARAM, "SerialParam")) {
              if (!is(BPPARAM, "SnowParam"))
                stop("At the moment parallelization only works with 'SnowParam' parallel backends.")

              if (BPPARAM$.clusterargs$type != "MPI")
                stop("At the moment the only type of parallel cluster available is 'MPI'") 

              if (BPPARAM$.clusterargs$spec > 0)
                clusterSize <- BPPARAM$.clusterargs$spec
            }

            rhs <- .parseFormula(model)

            pvaluesG0 <- nrr <- as(Matrix(numeric(), nrow=0, ncol=0), "dspMatrix")
            fix.Q <- restrict.Q <- NULL
            qorders <- integer()
            if (length(rhs$conditioning) > 0) {
              maskq0 <- sapply(rhs$conditioning, function(x) all(x$q==0))
              fix.Q <- unlist(lapply(rhs$conditioning[maskq0], '[[', "v"))
              restrict.Q <- unlist(lapply(rhs$conditioning[!maskq0], '[[', "v"))
              qorders <- lapply(rhs$conditioning[!maskq0], function(x) x$q)

              ## check that q-orders associated to different terms are identical
              if (!all(sapply(qorders, function(x, r) identical(x, r), qorders[[1]])))
                stop("q values should be identical throughout the given conditioning terms")

              qorders <- sort(qorders[1][[1]]) ## the [1][[1]] notation handles list() and [[1]] alone does not

              if (length(intersect(qorders, estimate@qOrders)) > 0)
                stop("The provided eQTLnetwork was estimated using already some of the given q-orders\n")

              if (is.null(qorders))
                qorders <- 0L

              if (any(qorders == 0)) {
                if (nrow(pvaluesG0) > 0)
                  stop("The provided eQTLnetwork has already a G^0 estimate\n")
                message("Calculating pairwise (conditional) independence tests\n")
              }
            }

            if (!identical(rhs$explanatory, unlist(lapply(.expandTerms(estimate@rhs$explanatory, param), function(x) x$v))))
              stop("Explanatory terms in the model formula do not match those of the given eQTLnetwork estimate\n")

            if (nrow(estimate@pvaluesG0 > 0)) #### CHECK THIS 
              pvaluesG0 <- estimate@pvaluesG0

            if (nrow(estimate@nrr > 0)) #### CHECK THIS
              nrr <- estimate@nrr

            res <- NULL
            if (length(rhs$conditioning) == 0 || any(qorders == 0)) {
              if (nrow(pvaluesG0) > 0)
                stop("The provided eQTLnetwork has already a G^0 estimate\n")
              rhs$conditioning <- estimate@rhs@conditioning

              message("Calculating pairwise (conditional) independence tests\n")
              pvaluesG0<- qpAllCItests(ggData(param), I=param@dVars, Q=NULL,
                                       pairup.i=geneNames(param), pairup.j=rhs$explanatory,
                                       return.type="p.value", clusterSize=clusterSize,
                                       verbose=verbose)$p.value
              qorders <- qorders[qorders > 0]
            }
           
            if (length(rhs$conditioning) > 0) {
              k <- length(estimate@qOrders)

              for (i in seq(along=qorders)) {
                q <- qorders[i]
                message(sprintf("Estimating non-rejection rates with q=%d\n", q))
                nrr.q <- qpNrr(ggData(param), I=param@dVars, q=q, pairup.i=geneNames(param), pairup.j=rhs$explanatory,
                               restrict.Q=restrict.Q, fix.Q=fix.Q, clusterSize=clusterSize, verbose=verbose)
                if (i+k == 1)
                  nrr <- nrr.q
                else
                  nrr <- (nrr.q + nrr*(i+k-1)) / (i+k)
              }
              qorders <- sort(c(qorders, estimate@qOrders))
              model <- .replaceFormulaQs(model, qorders)
            }

            new("eQTLnetwork", geneticMap=geneticMap(param), physicalMap=physicalMap(param),
                organism=param@organism, genome=param@genome, geneAnnotation=geneAnnotation(param),
                geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
                pvaluesG0=pvaluesG0, nrr=nrr, modelFormula=model, rhs=rhs, qOrders=qorders)
          })

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

.parseFormula <- function(f) {
  tryCatch({
    rhs <- deparse(f[[2]])
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

  rhs$explanatory <- unlist(lapply(.expandTerms(rhs$explanatory, param), function(x) x$v))
  if (length(rhs$conditioning) > 0)
    rhs$conditioning <- .expandTerms(rhs$conditioning, param)

  rhs
}

.replaceFormulaQs <- function(f, qorders) {
  tryCatch({
    rhs <- deparse(f[[2]])
  }, error=function(err) {
    cat("Malformed model formula\n")
    stop(conditionMessage(err), call.=FALSE)
  })

  rhs <- strsplit(rhs, " *\\| *")[[1]]
  if (length(rhs) > 2)
    stop("Conditioning operator '|' can only occur once in the model formula.")

  if (length(rhs) < 2)
    return(f)

  rhs <- strsplit(rhs, " *\\+ *")
  if (length(rhs) == 1)
    return(f)

  names(rhs) <- c("explanatory" ,"conditioning")

  rhs$conditioning <- sapply(rhs$conditioning,
                             function(term, qorders) {
                               m <- regexec("([^\\(]+)\\(q = ([0-9, ]+)\\)", term)
                               m <- regmatches(term, m)[[1]]
                               if (length(m) > 0) {
                                 term <- sprintf("%s(q = %s)", m[2], paste(qorders, collapse=", "))
                               }
                               term
                             }, qorders, USE.NAMES=FALSE)

  formula(paste("~", paste(paste(rhs$explanatory, collapse="+"),
                           paste(rhs$conditioning, collapse="+"),
                           sep="|")))
}
