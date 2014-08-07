## eQTLnetwork show method
setMethod("show", signature(object="eQTLnetwork"),
          function(object) {
            cat("eQTLnetwork object:\n")
            if (length(object@organism) > 0)
              cat(sprintf("  Organism: %s\n", object@organism))
            if (length(object@genome) > 0)
              cat(sprintf("  Genome: %s\n", object@genome))
            if (length(object@geneAnnotationTable) > 0)
              cat(sprintf("  Gene annotation: %s\n", object@geneAnnotationTable))
            mNames <- NULL
            if (length(object@geneticMap) > 0)
              mNames <- unlist(sapply(object@geneticMap, names), use.names=FALSE)
            else if (length(object@physicalMap) > 0)
              mNames <- unlist(sapply(object@physicalMap, names), use.names=FALSE)
            else
              warning("Both the genetic and the physical map are empty.")

            gNames <- names(object@geneAnnotation)
            cat(sprintf("  Input size: %d markers %d genes\n", length(mNames), length(gNames)))
            cat("  Model formula: ", deparse(object@modelFormula), "\n")
            g <- object@qpg@g
            if (numNodes(g) > 0) {
              edg <- extractFromTo(g)
              qstr <- NULL
              if (!is.na(object@p.value))
                qstr <- "0"
              if (!is.na(object@epsilon))
                qstr <- ifelse(!is.null(qstr), ## must be "0"
                               sprintf("0,%s", paste(object@qOrders, collapse=",")),
                               paste(object@qOrders, collapse=","))
              if (!is.na(object@alpha))
                qstr <- paste(qstr, "*", sep=",")
              padstr <- paste(rep(" ", times=nchar(qstr)+8), collapse="")
              cat(sprintf("  G^(%s): %d vertices and %d edges corresponding to\n",
                          qstr, numNodes(g), numEdges(g)))
              cat(sprintf("%s%d eQTL and %d gene-gene associations",
                          padstr, sum(edg$from %in% mNames | edg$to %in% mNames),
                          sum(edg$from %in% gNames & edg$to %in% gNames)))
              commaflag <- FALSE
              if (!is.na(object@p.value)) {
                cat(sprintf("meeting\n%sa %s-adjusted p-value < %.2f",
                            padstr, object@adjustMethod, object@p.value))
                commaflag <- TRUE
              }

              if (!is.na(object@epsilon)) {
                if (commaflag) cat(",\n") else cat("meeting\n")
                cat(sprintf("%sa non-rejection rate epsilon < %.2f",
                            padstr, object@epsilon))
                commaflag <- TRUE
              }

              if (!is.na(object@alpha)) {
                if (commaflag) cat(",\n") else cat("meeting\n")
                cat(sprintf("%sa forward eQTL selection significance level alpha < %.2f",
                            padstr, object@alpha))
              }

              cat(sprintf("\n%sand involving %d genes and %d eQTLs\n",
                          padstr, length(intersect(nodes(g), gNames)),
                          length(intersect(nodes(g), mNames))))
            }
          })

## alleQTL method
setMethod("alleQTL", signature(x="eQTLnetwork"),
          function(x) {
            mNames <- mChr <- mLoc <- NULL
            if (length(x@geneticMap) > 0) {
              mNames <- unlist(sapply(x@geneticMap, names), use.names=FALSE)
              mLoc <- unlist(x@geneticMap, use.names=FALSE)
              names(mLoc) <- mNames
              elen <- sapply(x@geneticMap, length)
              mChr <- rep(names(x@geneticMap), times=elen)
              names(mChr) <- mNames
            } else if (length(x@physicalMap) > 0) {
              mNames <- unlist(sapply(x@physicalMap, names), use.names=FALSE)
              mLoc <- unlist(x@physicalMap, use.names=FALSE)
              names(mLoc) <- mNames
              elen <- sapply(x@physicalMap, length)
              mChr <- rep(names(x@geneticMap), times=elen)
              names(mChr) <- mNames
            } else
              stop("Both the genetic and the physical map are empty.")

            gNames <- names(x@geneAnnotation)
            eQTLedges <- edges(x@qpg@g)[mNames]

            alledg <- extractFromTo(x@qpg@g)
            alledg$from <- as.character(alledg$from)
            alledg$to <- as.character(alledg$to)
            ggedg <- alledg[alledg$from %in% gNames & alledg$to %in% gNames, ]
            mgedg <- alledg[alledg$from %in% mNames | alledg$to %in% mNames, ]
            maskSwappedGenesMarkers <- mgedg$from %in% gNames
            if (any(maskSwappedGenesMarkers)) { ## put markers in 'from' and genes in 'to
              swappedGeneNames <- mgedg$from[maskSwappedGenesMarkers]
              mgedg$from[swappedGenesMarkers] <- mgedg$to[swappedGenesMarkers]
              mgedg$to[swappedGenesMarkers] <- mgedg$to[swappedGenesMarkers]
            }
            df <- as.data.frame(matrix(NA, nrow=0, ncol=4, dimnames=list(NULL, c("chrom", "location", "QTL", "gene"))))
            n.eqtl <- nrow(mgedg)
            if (n.eqtl > 0) {
              df <- data.frame(chrom=rep(NA_character_, n.eqtl),
                               location=rep(NA_real_, n.eqtl),
                               QTL=mgedg$from,
                               gene=mgedg$to,
                               stringsAsFactors=FALSE)
              df$chrom <- mChr[df$QTL]
              df$location <- mLoc[df$QTL]
            }
            df
          })

