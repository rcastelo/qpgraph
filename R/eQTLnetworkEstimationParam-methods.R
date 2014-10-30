.isPkgLoaded <- function(name) {
  (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}

## eQTLnetworkEstimationParam constructor
eQTLnetworkEstimationParam <- function(ggData, geneticMap=NULL, physicalMap=NULL,
                                       dVars=NULL, genes=NULL, geneAnnotation, genome=Seqinfo()) {

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
    markerNames <- unlist(lapply(geneticMap, names), use.names=FALSE)
  } else {                                             ## input data is a 'matrix' or 'data.frame' object
    if (!is.null(dVars)) {
      if (any(!dVars %in% colnames(ggData)))
        stop("Some variables in 'dVars' do not form part of the input data in 'ggData'.")
    }

    if (!is.null(geneticMap)) {
      if (!is(geneticMap, "list") && !is(geneticMap, "map"))
        stop("'geneticMap' should be either a 'qtl::map' object or a 'list' object.")
      markerNames <- unlist(lapply(geneticMap, names), use.names=FALSE)
    }

    if (!is.null(physicalMap)) {
      if (!is(physicalMap, "list") && !is(physicalMap, "map"))
        stop("'physicalMap' should be either a 'qtl::map' object or a 'list' object.")
      markerNames <- unlist(lapply(physicalMap, names), use.names=FALSE)
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

  if (missing(geneAnnotation))
    stop("argument 'geneAnnotation' must be set either as a character string, a 'data.frame' object or a TxDb object.")

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
    if (!is(geneAnnotation, "TxDb"))
      stop(sprintf("the object loaded with name %s is not a 'TxDb' object.", geneAnnotation))
  } else if (class(geneAnnotation) == "data.frame") {
    if (is.null(rownames(geneAnnotation)))
      stop("when argument 'geneAnnotation' is a 'data.frame' object, it should have row names uniquely identifying each gene.")
    if (!is(genome, "Seqinfo"))
      stop("when argument 'geneAnnotation' is a 'data.frame' object, argument 'genome' should contain a 'Seqinfo' object.")

    clsannot <- sapply(geneAnnotation, class)
       
    if ((clsannot[1] != "character" && clsannot[1] != "factor") ||
        (clsannot[2] != "numeric" && clsannot[2] != "integer") ||
        (clsannot[3] != "numeric" && clsannot[3] != "integer") ||
        (clsannot[4] != "character" && clsannot[4] != "factor"))
      stop("the first four column classes in the 'data.frame' object 'geneAnnotation' must be character or factor, numeric or integer, numeric or integer and character or factor, corresponding to chromosome, start position, end position and strand of the gene annotation in the genome.")
    geneAnnotation[[1]] <- as.character(geneAnnotation[[1]])
    geneAnnotation[[4]] <- as.character(geneAnnotation[[4]])

    map <- geneticMap
    if (is.null(map))
      map <- physicalMap

    geneIDs <- rownames(geneAnnotation)

    geneAnnotation <- GRanges(seqnames=geneAnnotation[[1]],
                              IRanges(start=geneAnnotation[[2]], end=geneAnnotation[[3]]),
                              strand=geneAnnotation[[4]],
                              seqinfo=genome)
    names(geneAnnotation) <- geneIDs
  } else if (!is(geneAnnotation, "TxDb"))
    stop("argument 'geneAnnotation' must be either a character string, a 'data.frame' object or a 'TxDb' object")

  organism <- geneAnnotationTable <- character()
  if (is(geneAnnotation, "TxDb")) {
    md <- metadata(geneAnnotation)
    organism <- md[md$name %in% "Organism", "value"]
    ## genome <- md[md$name %in% "Genome", "value"]
    genome <- seqinfo(geneAnnotation)
    geneAnnotationTable <- paste(md[grep("Table", md$name)[1], ], collapse=" ") 

    geneAnnotation <- transcripts(geneAnnotation)
    names(geneAnnotation) <- geneAnnotation$tx_name
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
  if (any(!genes %in% names(geneAnnotation))) {
    maskmissinggenes <- !genes %in% names(geneAnnotation)
    if (any(genes[maskmissinggenes] %in% make.names(names(geneAnnotation))))
      names(geneAnnotation) <- make.names(names(geneAnnotation))

    if (any(!genes %in% names(geneAnnotation))) {
      warning(sprintf("%d genes were not found in 'geneAnnotation'. They will be removed from the input data for further analysis.",
                      sum(!genes %in% names(geneAnnotation))))
      mt <- match(genes[!genes %in% names(geneAnnotation)], colnames(ggDataMatrix))
      ggDataMatrix <- ggDataMatrix[, -mt, drop=FALSE]
      genes <- genes[genes %in% names(geneAnnotation)]
    }
  }

  ## remove annotations from genes not in the data
  geneAnnotation <- geneAnnotation[genes]

  
  new("eQTLnetworkEstimationParam", ggData=ggDataMatrix, geneticMap=geneticMap,
      physicalMap=physicalMap, organism=organism, genome=genome,
      geneAnnotation=geneAnnotation, geneAnnotationTable=geneAnnotationTable,
      dVars=dVars)
}

## eQTLnetworkEstimationParam show method
setMethod("show", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            cat("eQTLnetworkEstimationParam object:\n")
            if (length(object@organism) > 0)
              cat(sprintf("  Organism: %s\n", object@organism))
            if (length(object@genome) > 0)
              cat(sprintf("  Genome: %s\n", unique(genome(object@genome))))
            if (length(object@geneAnnotationTable) > 0)
              cat(sprintf("  Gene annotation: %s\n", object@geneAnnotationTable))
            cat(sprintf("  %d markers %d genes\n",
                max(c(length(unlist(object@geneticMap, use.names=FALSE)),
                      length(unlist(object@physicalMap, use.names=FALSE)))),
                length(object@geneAnnotation)))

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

setMethod("ggData", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@ggData
          })

setMethod("markerNames", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            markers <- character()
            if (!is.null(object@geneticMap))
              markers <- unlist(lapply(object@geneticMap, names), use.names=FALSE)

            if (!is.null(object@physicalMap))
              markers <- unlist(lapply(object@physicalMap, names), use.names=FALSE)

            markers
          })

setMethod("geneNames", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            names(object@geneAnnotation)
          })

setMethod("geneAnnotation", signature(object="eQTLnetworkEstimationParam"),
          function(object) {
            object@geneAnnotation
          })
