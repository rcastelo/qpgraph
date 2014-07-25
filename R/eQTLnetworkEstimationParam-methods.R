.isPkgLoaded <- function(name) {
  (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}

## eQTLnetworkEstimationParam constructor
eQTLnetworkEstimationParam <- function(ggData, geneticMap=NULL, physicalMap=NULL,
                                       dVars=NULL, genes=NULL, geneAnnotation,
                                       qOrders=NA_integer_) {

  ## check that 'ggData' is one of the allowed object classes
  if (!is(ggData, "cross") && !is(ggData, "matrix") && is(ggData, "data.frame"))
    stop("Argument 'ggData' must be a either a 'qtl::cross' object, a 'matrix' object or a 'data.frame' object.")

  ## if 'ggData' is a cross object, then there is no need for the 'geneticMap' argument
  if (is(ggData, "cross") && !is.null(geneticMap))
    stop("If argument 'ggData' is a 'qtl::cross' object, then argument 'geneticMap' must not be set.")

  ## if 'ggData' is a 'matrix' or a 'data.frame', then we need either a 'geneticMap' or a 'physicalMap'
  if ((is(ggData, "matrix") || is(ggData, "data.frame")) &&
      is.null(geneticMap) && is.null(physicalMap))
    stop("When 'ggData' is a 'matrix' or a 'data.frame', then either 'geneticMap' or 'physicalMap' must be set.")

  markerNames <- NULL
  if (is(ggData, 'cross')) {                           ## input data is a 'cross' object
    if (!is.null(dVars))
      stop("When 'ggData' is a 'qtl::cross' object, then argument 'dVars' must not be set.")

    geneticMap <- pull.map(ggData)
    markerNames <- names(unlist(geneticMap))
  } else {                                             ## input data is a 'matrix' or 'data.frame' object
    if (!is.null(dVars)) {
      if (any(!dVars %in% colnames(ggData)))
        stop("Some variables in 'dVars' do not form part of the input data in 'ggData'.")
    }

    if (!is.null(geneticMap)) {
      if (!is(geneticMap, "list") && !is(geneticMap, "map"))
        stop("'geneticMap' should be either a 'qtl::map' object or a 'list' object.")
      markerNames <- names(unlist(geneticMap))
    }

    if (!is.null(physicalMap)) {
      if (!is(physicalMap, "list") && !is(physicalMap, "map"))
        stop("'physicalMap' should be either a 'qtl::map' object or a 'list' object.")
      markerNames <- names(unlist(physicalMap))
    }

    if (any(!markerNames %in% colnames(ggData)))
      stop("Some markers do not form part of the input data in 'ggData'.")
  }

  if (length(markerNames) == 0)
    stop("no markers present in the input data in 'ggData'.")

  if (is.null(genes)) { ## assume everything except markers and discrete variables are genes

    if (is(ggData, "cross")) {
      dfclass <- sapply(ggData$pheno, class)
      discreteMask <- dfclass == "character" | dfclass == "factor" | dfclass == "logical"
      genes <- setdiff(colnames(ggData$pheno), c(markerNames, colnames(ggData$pheno)[discreteMask]))
    }

    if (is(ggData, "matrix")) {
      if (any(!dVars %in% colnames(ggData)))
        stop("Some of the discrete variables specified in 'dVars' do not exist in the input data 'ggData'.")
      genes <- setdiff(colnames(ggData), c(markerNames, dVars))
    }

    if (is(ggData, "data.frame")) {
      dfclass <- sapply(ggData, class)
      discreteMask <- dfclass == "character" | dfclass == "factor" | dfclass == "logical"
      genes <- setdiff(colnames(ggData), c(markerNames, colnames(ggData)[discreteMask]))
    }

  } else { ## genes are specified

    if (is(ggData, "cross")) {
      if (any(!genes %in% colnames(ggData$pheno)))
        stop("Some identifiers in argument 'genes' do not form part of the input data in 'ggData'.")
    } else {
      if (any(!genes %in% colnames(ggData)))
        stop("Some identifiers in argument 'genes' do not form part of the input data in 'ggData'.")
    }

  }

  if (length(genes) == 0)
    stop("no genes present in the input data in 'ggData'.")

  if (class(geneAnnotation) == "character") {
    if (!exists(geneAnnotation)) {
      if (!geneAnnotation %in% installed.packages()[, "Package"])
        stop(sprintf("Please install the Bioconductor package %s.", geneAnnotation))

      if (!.isPkgLoaded(geneAnnotation)) {
        message("Loading TxDb annotation package ", geneAnnotation)
        if (!suppressPackageStartupMessages(require(geneAnnotation, character.only=TRUE)))
          stop(sprintf("package %s could not be loaded.", geneAnnotation))
      }
    }
    tryCatch({
      geneAnnotation <- get(geneAnnotation)
    }, error=function(err) {
      stop(sprintf("package %s does not load a object with the same name of the package."))
    })
    if (!is(geneAnnotation, "TranscriptDb"))
      stop(sprintf("The object loaded with name %s is not a 'TranscriptDb' object.", geneAnnotation))
  }

  organism <- genome <- geneAnnotationTable <- character()
  if (is(geneAnnotation, "TranscriptDb")) {
    md <- metadata(geneAnnotation)
    organism <- md[md$name %in% "Organism", "value"]
    genome <- md[md$name %in% "Genome", "value"]
    geneAnnotationTable <- paste(md[grep("Table", md$name)[1], ], collapse=" ") 

    chr2idx <- seqlevels(geneAnnotation)
    chr2idx <- do.call("names<-", list(orderSeqlevels(chr2idx), chr2idx))
    suppressWarnings(geneAnnotation <- select(geneAnnotation, keys=keys(geneAnnotation),
                                              columns=c("GENEID", "TXCHROM", "TXSTART"),
                                              keytype="GENEID"))
    geneAnnotation$GENEID <- make.names(geneAnnotation$GENEID)
    minStart <- split(geneAnnotation$TXSTART, geneAnnotation$GENEID)
    minStart <- sapply(minStart, min)
    chr <- split(geneAnnotation$TXCHROM, geneAnnotation$GENEID)
    chr <- sapply(chr, unique)
    chr <- do.call("names<-", list(chr2idx[chr], names(chr)))
    geneAnnotation <- cbind(as.integer(chr), as.integer(minStart))
    rownames(geneAnnotation) <- names(chr)
  }

  ## store the all the input data into a big matrix
  ## this should be in the near future efficiently handled
  ## to reduce the memory footprint
  ggDataMatrix <- ggData
  if (is(ggData, "cross")) {
    phenoNames <- colnames(ggData$pheno)
    phclass <- sapply(ggData$pheno, class)
    ggDataMatrix <- matrix(NA_real_, nrow=qtl::nind(ggData),
                           ncol=qtl::nphe(ggData)+qtl::totmar(ggData),
                           dimnames=list(NULL, c(markerNames, phenoNames)))
    ggDataMatrix[, markerNames] <- do.call("cbind", lapply(ggData$geno, function(x) x$data))
    if (any(phclass == "numeric"))
      ggDataMatrix[, phenoNames[phclass == "numeric"]] <-
        as.matrix(ggData$pheno[, phenoNames[phclass == "numeric"]])

    if (any(phclass == "integer"))
      ggDataMatrix[, phenoNames[phclass == "integer"]] <-
        do.call("mode<-", list(as.matrix(ggData$pheno[, phenoNames[phclass == "integer"]]), "numeric"))
    
    discreteMask <- phclass == "character" | phclass == "factor" | phclass == "logical"
    dVars <- c(markerNames, phenoNames[discreteMask])
    if (any(discreteMask))
      ggDataMatrix[, phenoNames[discreteMask]] <-
        as.matrix(as.data.frame(lapply(lapply(ggData$pheno[, phenoNames[discreteMask]], factor), as.numeric)))
  } else if (is(ggData, "data.frame")) {
    cnames <- colnames(ggData)
    cclass <- sapply(ggData, class)
    ggDataMatrix <- matrix(NA_real_, nrow=nrow(ggData), ncol=ncol(ggData),
                           dimnames=list(NULL, c(markerNames, setdiff(cnames, markerNames))))
    discreteMask <- cclass == "character" | cclass == "factor" | cclass == "logical"
    dVars <- unique(c(markerNames, cnames[discreteMask]))
    ggDataMatrix[, dVars] <- 
        as.matrix(as.data.frame(lapply(lapply(ggData[, dVars], factor), as.numeric)))

    if (any(cclass == "numeric"))
      ggDataMatrix[, cnames[cclass == "numeric"]] <-
        as.matrix(ggData[, cnames[cclass == "numeric"]])

    if (any(cclass == "integer"))
      ggDataMatrix[, cnames[cclass == "integer"]] <-
        do.call("mode<-", list(as.matrix(ggData[, cnames[cclass == "integer"]]), "numeric"))
  }

  ## remove from data genes without annotation
  if (any(!genes %in% rownames(geneAnnotation))) {
    warning(sprintf("%d genes were not found in 'geneAnnotation'. They will be removed from the input data for further analysis.",
                    sum(!genes %in% rownames(geneAnnotation))))
    mt <- match(genes[!genes %in% rownames(geneAnnotation)], colnames(ggDataMatrix))
    ggDatamatrix <- ggDataMatrix[, -mt]
    genes <- genes[genes %in% rownames(geneAnnotation)]
    geneAnnotation <- geneAnnotation[genes, ]
  }

  
  new("eQTLnetworkEstimationParam", ggData=ggDataMatrix, geneticMap=geneticMap,
      physicalMap=physicalMap, organism=organism, genome=genome,
      geneAnnotation=geneAnnotation, geneAnnotationTable=geneAnnotationTable,
      dVars=dVars, qOrders=as.integer(qOrders))
}

## eQTLnetworkEstimationParam show method
setMethod("show", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            cat("eQTLnetworkEstimationParam object:\n")
            if (length(object@organism) > 0)
              cat(sprintf("  Organism: %s\n", object@organism))
            if (length(object@genome) > 0)
              cat(sprintf("  Genome: %s\n", object@genome))
            if (length(object@geneAnnotationTable) > 0)
              cat(sprintf("  Gene annotation: %s\n", object@geneAnnotationTable))
            cat(sprintf("  %d markers %d genes\n",
                max(c(length(unlist(object@geneticMap, use.names=FALSE)),
                      length(unlist(object@physicalMap, use.names=FALSE)))),
                nrow(object@geneAnnotation)))
            cat(sprintf("  q orders: {%s}\n", paste(object@qOrders, collapse=", ")))

            invisible(object)
          })

## setter and getter methods
setMethod("geneticMap", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@geneticMap
          })

setMethod("physicalMap", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@physicalMap
          })

setMethod("geneAnnotation", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@physicalMap
          })

setMethod("qOrders", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@qOrders
          })

setReplaceMethod("qOrders", signature(object="eQTLnetworkEstimationParam", value="ANY"),
                 function(object) {
                   object@qOrders <- as.integer(value)
                 })
