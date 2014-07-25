.isPkgLoaded <- function(name) {
  (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}

## eQTLnetworkEstimationParam constructor
eQTLnetworkEstimationParam <- function(ggData, geneticMap=NULL, physicalMap=NULL, genes=NULL, geneAnnotation,
                                       nullHypothesis=c("noeQTLanywhere", "noeQTLatMarker"),
                                       pValueCutoff=0.05, adjustMethod=p.adjust.methods, qOrders=NA_integer_) {
  nullHypothesis <- match.arg(nullHypothesis)

  ## check that 'ggData' is one of the allowed object classes
  if (!is(ggData, "cross") && !is(ggData, "matrix") && is(ggData, "data.frame"))
    stop("Argument 'ggData' must be a either a 'qtl::cross' object, a 'matrix' object or a 'data.frame' object.")

  ## if 'ggData' is a cross object, then there is no need for the 'geneticMap' argument
  if (is(ggData, "cross") && !is.null(geneticMap))
    stop("If argument 'ggData' is a 'cross' object, then argument 'geneticMap' must not be set.")

  ## if 'ggData' is a 'matrix' or a 'data.frame', then we need either a 'geneticMap' or a 'physicalMap'
  if ((is(ggData, "matrix") || is(ggData, "data.frame")) &&
      is.null(geneticMap) && is.null(physicalMap))
    stop("When 'ggData' is a 'matrix' or a 'data.frame', then either 'geneticMap' or 'physicalMap' must be set.")

  markerNames <- NULL
  if (is(ggData, 'cross')) {                           ## input data is a 'cross' object
    geneticMap <- pull.map(ggData)
    markerNames <- names(unlist(geneticMap))
  } else {                                             ## input data is a 'matrix' or 'data.frame' object
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
      stop("Some markers do not form part of the input data.")
  }

  if (length(markerNames) == 0)
    stop("no markers present in the input data.")

  if (is.null(genes)) {
    if (is(ggData, "cross"))
      genes <- colnames(ggData$pheno)

    if (is(ggData, "matrix") || is(ggData, "data.frame"))
      genes <- setdiff(colnames(ggData), markerNames)
  } else {
    if (is(ggData, "cross"))
      if (any(!genes %in% colnames(ggData$pheno)))
        stop("Some identifiers in argument 'genes' do not form part of the input data in 'ggData'.")

      if (any(!genes %in% colnames(ggData)))
        stop("Some identifiers in argument 'genes' do not form part of the input data in 'ggData'.")
  }

  if (length(genes) == 0)
    stop("no markers present in the input data.")

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

  if (is(geneAnnotation, "TranscriptDb")) {
      geneAnnotation <- select(geneAnnotation, keys=keys(geneAnnotation),
                               columns=c("GENEID", "TXCHROM", "TXSTART"),
                               keytype="GENEID")
      minStart <- split(geneAnnotation$TXSTART, geneAnnotation$GENEID)
      minStart <- sapply(minStart, min)
      chr <- split(geneAnnotation$TXCHROM, geneAnnotation$GENEID)
      chr <- sapply(chr, unique)
      chr <- gsub("chr", "", chr)
      chr <- gsub("M", "XVII", chr)
      chr <- as.roman(chr)
      geneAnnotation <- cbind(as.integer(chr), as.integer(minStart))
      rownames(geneAnnotation) <- names(chr)
  }
  
  new("eQTLnetworkEstimationParam", geneticMap=geneticMap, physicalMap=physicalMap,
      geneAnnotation=geneAnnotation, nullHypothesis=nullHypothesis, adjustMethod=adjustMethod,
      pValueCutoff=pValueCutoff, qOrders=qOrders)
}

## eQTLnetworkEstimationParam show method
setMethod("show", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            cat("eQTLnetworkEstimationParam object:\n")
            cat(sprintf("  %d markers %d genes\n",
                max(c(length(unlist(object@geneticMap, use.names=FALSE)),
                      length(unlist(object@physicalMap, use.names=FALSE)))),
                nrow(object@geneAnnotation)))
            cat(sprintf("  Null hypothesis for G^(0): %s\n", object@nullHypothesis))
            if (object@nullHypothesis == "noeQTLatMarker")
              cat(sprintf("  P-value adjusting method for G^(0): %s\n", object@adjustMethod))
            cat(sprintf("  P-value cutoff for G^(0): %.2f\n", object@pValueCutoff))
            cat(sprintf("  q orders: {%s}\n", paste(object@qOrders, collapse=", ")))

            invisible(object)
          })
