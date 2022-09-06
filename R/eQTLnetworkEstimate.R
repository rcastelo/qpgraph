## methods for estimating eQTL networks

## estimate parameters of the eQTLnetwork
setMethod("eQTLnetworkEstimate", signature=c(param="eQTLnetworkEstimationParam",
                                             model="formula",
                                             estimate="missing"),
          function(param, model, estimate, p.value=NA_real_, method=p.adjust.methods,
                   epsilon=NA_real_, alpha=NA_real_, verbose=TRUE,
                   BPPARAM=bpparam("SerialParam")) {

            if (!is.na(p.value) || !is.na(epsilon))
              stop("'p.value' or 'epsilon' cutoffs can only be specified without 'model'.")

            clusterSize <- 1 ## this should change to use BiocParallel
            if (!is(BPPARAM, "SerialParam")) {
              if (!is(BPPARAM, "SnowParam"))
                stop("At the moment parallelization only works with 'SnowParam' parallel backends.")

              if (BPPARAM$.clusterargs$type != "MPI")
                stop("At the moment the only type of parallel cluster available is 'MPI'") 

              if (BPPARAM$.clusterargs$spec > 0)
                clusterSize <- BPPARAM$.clusterargs$spec
            }

            rhs <- .parseFormula(model, param)

            pvaluesG0 <- nrr <- new("dspMatrix", x=numeric())
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
                                        long.dim.are.variables=FALSE, exact.test=TRUE,
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
                               long.dim.are.variables=FALSE, exact.test=TRUE,
                               restrict.Q=restrict.Q, fix.Q=fix.Q,
                               clusterSize=clusterSize, verbose=verbose)
                if (i+k == 1)
                  nrr <- nrr.q
                else
                  nrr <- (nrr.q + nrr*(i+k-1)) / (i+k)
              }
            }

            g.0 <- graphBAM(df=data.frame(from=character(), to=character(), weight=integer()))
            qpg <- new("qpGraph", p=0L, q=integer(), n=NA_integer_, epsilon=NA_real_, g=g.0)
            new("eQTLnetwork", geneticMap=geneticMap(param),
                physicalMap=physicalMap(param), organism=param@organism,
                genome=param@genome, geneAnnotation=geneAnnotation(param),
                geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
                pvaluesG0=pvaluesG0, nrr=nrr, modelFormula=model, rhs=rhs,
                qOrders=qorders, p.value=NA_real_, adjustMethod="none",
                epsilon=NA_real_, alpha=NA_real_, qpg=qpg)
          })

## estimate parameters of the eQTLnetwork given an existing estimate
setMethod("eQTLnetworkEstimate", signature=c(param="eQTLnetworkEstimationParam",
                                             model="formula",
                                             estimate="eQTLnetwork"),
          function(param, model, estimate, p.value=NA_real_, method=p.adjust.methods,
                   epsilon=NA_real_, alpha=NA_real_, verbose=TRUE,
                   BPPARAM=bpparam("SerialParam")) {

            if (!is.na(p.value) || !is.na(epsilon) || !is.na(alpha))
              stop("'p.value', 'epsilon' or 'alpha' cutoffs can only be specified without 'model'.")

            clusterSize <- 1 ## this should change to use BiocParallel
            if (!is(BPPARAM, "SerialParam")) {
              if (!is(BPPARAM, "SnowParam"))
                stop("At the moment parallelization only works with 'SnowParam' parallel backends.")

              if (BPPARAM$.clusterargs$type != "MPI")
                stop("At the moment the only type of parallel cluster available is 'MPI'") 

              if (BPPARAM$.clusterargs$spec > 0)
                clusterSize <- BPPARAM$.clusterargs$spec
            }

            rhs <- .parseFormula(model, param)

            pvaluesG0 <- nrr <- new("dspMatrix", x=numeric())
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

            ## take p-value, adjust method and qp-graph parameters from the estimate
            p.value <- estimate@p.value
            method <- estimate@adjustMethod
            qpg <- estimate@qpg

            if (nrow(estimate@nrr > 0)) #### CHECK THIS
              nrr <- estimate@nrr

            res <- NULL
            if (length(rhs$conditioning) == 0 || any(qorders == 0)) {
              if (nrow(pvaluesG0) > 0)
                stop("The provided eQTLnetwork has already a G^0 estimate\n")
              rhs$conditioning <- estimate@rhs@conditioning

              message("Calculating pairwise (conditional) independence tests\n")
              pvaluesG0 <- qpAllCItests(ggData(param), I=param@dVars, Q=NULL,
                                        pairup.i=geneNames(param), pairup.j=rhs$explanatory,
                                        long.dim.are.variables=FALSE, exact.test=TRUE,
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
                               long.dim.are.variables=FALSE, exact.test=TRUE,
                               restrict.Q=restrict.Q, fix.Q=fix.Q,
                               clusterSize=clusterSize, verbose=verbose)
                if (i+k == 1)
                  nrr <- nrr.q
                else
                  nrr <- (nrr.q + nrr*(i+k-1)) / (i+k)
              }
              qorders <- sort(c(qorders, estimate@qOrders))
              model <- .replaceFormulaQs(model, qorders)
            }

            ## g.0 <- graphBAM(df=data.frame(from=character(), to=character(), weight=integer()))
            ## qpg <- new("qpGraph", p=0L, q=integer(), n=NA_integer_, epsilon=NA_real_, g=g.0)

            new("eQTLnetwork", geneticMap=geneticMap(param),
                physicalMap=physicalMap(param), organism=param@organism,
                genome=param@genome, geneAnnotation=geneAnnotation(param),
                geneAnnotationTable=param@geneAnnotationTable, dVars=param@dVars,
                pvaluesG0=pvaluesG0, nrr=nrr, modelFormula=model, rhs=rhs,
                qOrders=qorders, p.value=p.value, adjustMethod=method,
                epsilon=epsilon, alpha=alpha, qpg=qpg)
          })


## estimate the eQTLnetwork based on the estimated parameters and given cutoffs
setMethod("eQTLnetworkEstimate", signature=c(param="eQTLnetworkEstimationParam",
                                             model="missing",
                                             estimate="eQTLnetwork"),
          function(param, model, estimate, p.value=NA_real_, method=p.adjust.methods,
                   epsilon=NA_real_, alpha=NA_real_, verbose=TRUE,
                   BPPARAM=bpparam("SerialParam")) {

            method <- match.arg(method)

            if (!is.na(p.value) && nrow(estimate@pvaluesG0) < 1)
              stop("argument 'estimate' has no pairwise (conditional) independence tests.")

            if (!is.na(epsilon) && nrow(estimate@nrr) < 1)
              stop("argument 'estimate' has no non-rejection rates.")

            mNames <- markerNames(estimate)
            if (!identical(mNames, markerNames(param)))
              stop("Markers are different between 'param' and 'estimate'.")

            gNames <- names(estimate@geneAnnotation)
            if (!identical(gNames, geneNames(param)))
              stop("Genes are different between 'param' and 'estimate'.")

            qpg <- estimate@qpg

            if (!is.na(p.value)) {
              ## here we assume that the diagonal is set to NA
              p.adj <- estimate@pvaluesG0[upper.tri(estimate@pvaluesG0, diag=TRUE)]
              p.adj[!is.na(p.adj)] <- p.adjust(p.adj[!is.na(p.adj)], method=method)
              p.adj <- new("dspMatrix", Dim=dim(estimate@pvaluesG0),
                           Dimnames=dimnames(estimate@pvaluesG0), x=p.adj)
              g.0 <- qpAnyGraph(p.adj, threshold=p.value, remove="above",
                                decreasing=FALSE)
              if (qpg@p > 0) {
                if (nrow(ggData(param)) != qpg@n)
                  stop("The qp-graph in the given eQTL network was estimated from different data.")
                g <- graphIntersect(qpg@g, g.0)
                qpg <- new("qpGraph", p=numNodes(g), q=unique(sort(c(0L, qpg@q))),
                           n=qpg@n, epsilon=qpg@epsilon, g=g)
              } else
                qpg <- new("qpGraph", p=numNodes(g.0), q=0L, n=nrow(ggData(param)),
                           epsilon=NA_real_, g=g.0)
            }

            if (!is.na(epsilon)) {
              if (nrow(estimate@nrr) < 1)
                stop("non-rejection rates have not been calculated.")

              if (is.na(p.value)) {
                p.value <- estimate@p.value
                method <- estimate@adjustMethod
              }

              g.q <- qpGraph(nrrMatrix=estimate@nrr, epsilon=epsilon, q=estimate@qOrders,
                             n=nrow(ggData(param)))
              if (qpg@p > 0) {
                g <- graphIntersect(qpg@g, g.q@g)
                qpg <- new("qpGraph", p=numNodes(g), q=unique(sort(c(qpg@q, estimate@qOrders))),
                           n=qpg@n, epsilon=epsilon, g=g)
              } else
                qpg <- g.q
            }

            if (!is.na(alpha)) {
              epsilon <- estimate@epsilon
              if (is.na(p.value)) {
                p.value <- estimate@p.value
                method <- estimate@adjustMethod
              }

              if (numNodes(qpg@g) < 1 || nrow(estimate@nrr) < 1)
                stop("forward selection with an 'alpha' cutoff requires non-rejection rates and a given 'epsilon' cutoff")
              if (nrow(estimate@nrr) < 1)
                stop("non-rejection rates have not been calculated.")

              edg <- edges(qpg@g)[intersect(nodes(qpg@g), gNames)]
              edg <- lapply(edg, function(x, I) intersect(x, I), mNames)
              edg <- edg[sapply(edg, length) > 0]
              edg <- mapply(function(m, g) c(g, m), edg, as.list(names(edg)),
                            SIMPLIFY=FALSE) ## add the gene before markers
              edg <- bplapply(edg,
                              function(v, X, nrr, alpha) {
                                o <- order(nrr[cbind(rep(v[1], times=length(v)-1), v[-1])])
                                v[-1] <- v[-1][o]
                                dropMask <- rep(FALSE, times=length(v[-1]))
                                Q <- NULL
                                pv <- 0
                                i <- 1
                                while (i <= length(v[-1]) && pv <= alpha) { ## fwd selection
                                  pv <- qpCItest(X, i=v[1], j=v[i+1], Q=Q, I=v[-1], long.dim.are.variables=FALSE)$p.value
                                  dropMask[i] <- pv > alpha
                                  Q <- c(Q, v[i+1])
                                  i <- i + 1
                                }
                                if (i < length(v[-1]))
                                  dropMask[i:length(v[-1])] <- TRUE
                                ## NAs may occur when not sufficient data is available
                                v[-1][dropMask | is.na(dropMask)]
                              }, ggData(param), estimate@nrr, alpha=alpha, BPPARAM=BPPARAM)
              edg <- edg[sapply(edg, length) > 0]
              if (length(edg) > 0) {
                elen <- sapply(edg, length)
                qpg@g <- removeEdge(from=rep(names(edg), times=elen), to=unlist(edg, use.names=FALSE), qpg@g)
                ## alledg <- extractFromTo(qpg@g)
                ## alledg$from <- as.character(alledg$from)
                ## alledg$to <- as.character(alledg$to)
                ## ggedg <- alledg[alledg$from %in% gNames & alledg$to %in% gNames, ]
                ## mgedg <- alledg[alledg$from %in% mNames | alledg$to %in% mNames, ]
                ## maskSwappedGenesMarkers <- mgedg$from %in% gNames
                ## if (any(maskSwappedGenesMarkers)) { ## put markers in 'from' and genes in 'to
                ##   swappedGeneNames <- mgedg$from[maskSwappedGenesMarkers]
                ##   mgedg$from[maskSwappedGenesMarkers] <- mgedg$to[maskSwappedGenesMarkers]
                ##   mgedg$to[maskSwappedGenesMarkers] <- swappedGeneNames
                ## }
                ## elen <- sapply(edg, length)
                ## rmedg <- data.frame(from=unlist(edg, use.names=FALSE),
                ##                     to=rep(names(edg), times=elen))
                ## keepMask <- !paste(mgedg$from, mgedg$to, sep="__") %in% paste(rmedg$from, rmedg$to, sep="__")
                ## mgedg <- mgedg[keepMask, ]
                ## qpg@g <- graphBAM(rbind(ggedg, mgedg, stringsAsFactors=FALSE))
              }
            }

            new("eQTLnetwork", geneticMap=geneticMap(param),
                physicalMap=physicalMap(param),
                organism=param@organism, genome=param@genome,
                geneAnnotation=geneAnnotation(param),
                geneAnnotationTable=param@geneAnnotationTable,
                dVars=param@dVars, pvaluesG0=estimate@pvaluesG0,
                nrr=estimate@nrr, modelFormula=estimate@modelFormula,
                rhs=estimate@rhs, qOrders=estimate@qOrders,
                p.value=p.value, adjustMethod=method,
                epsilon=epsilon, alpha=alpha, qpg=qpg)
          })

setMethod("varExplained", signature=c(param="eQTLnetworkEstimationParam",
                                      estimate="eQTLnetwork"),
          function(param, estimate, BPPARAM=bpparam("SerialParam")) {
            if (nrow(estimate@nrr) < 1)
              stop("non-rejection rates are required at the moment to calculate the explained variance per gene\n")
            eqtls <- alleQTL(estimate)
            eqtls <- split(eqtls$QTL, eqtls$gene)
            ## add gene before markers
            eqtls <- mapply(function(m, g) c(g, m), eqtls, as.list(names(eqtls)))
            eta2xgene <- bplapply(eqtls,
                                  function(qtls, X, nrr) {
                                    o <- order(nrr[cbind(rep(qtls[1], times=length(qtls)-1), qtls[-1])])
                                    qtls[-1] <- qtls[-1][o]
                                    eta2 <- 0
                                    Q <- NULL
                                    for (i in 1:length(qtls[-1])) {
                                      this.eta2 <- qpCItest(X, i=qtls[1], j=qtls[i+1], Q=Q, I=qtls[-1], long.dim.are.variables=FALSE)$estimate
                                      eta2 <- eta2 + this.eta2
                                      Q <- c(Q, qtls[i+1])
                                    }
                                    eta2
                                  }, ggData(param), estimate@nrr, BPPARAM=BPPARAM)
            eta2xgene <- unlist(lapply(eta2xgene, as.vector))
            df <- data.frame(eta2eQTLs=eta2xgene)
            rownames(df) <- names(eta2xgene)
            df
          })

##
## private functions
##

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

.parseFormula <- function(f, param, expand=TRUE) {
  tryCatch({
    rhs <- paste(deparse(f[[2]]), collapse="")
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

  if (expand) {
    rhs$explanatory <- unlist(lapply(.expandTerms(rhs$explanatory, param), function(x) x$v))
    if (length(rhs$conditioning) > 0)
      rhs$conditioning <- .expandTerms(rhs$conditioning, param)
  }

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
