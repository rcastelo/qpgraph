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
                cat(sprintf(" meeting\n%sa %s-adjusted p-value < %.2f",
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

## getter and setter methods
setMethod("geneticMap", signature(object="eQTLnetwork"),
          function(object) {
            object@geneticMap
          })

setMethod("physicalMap", signature(object="eQTLnetwork"),
          function(object) {
            object@physicalMap
          })

setMethod("geneAnnotation", signature(object="eQTLnetwork"),
          function(object) {
            object@geneAnnotation
          })

setMethod("alleQTL", signature(x="eQTLnetwork"),
          function(x, map=c("genetic", "physical")) {
            map <- match.arg(map)
            gMap <- geneticMap(x)
            pMap <- physicalMap(x)

            if (length(gMap) < 1 && length(pMap) < 1)
              stop("Both, genetic and physical maps in the input object are empty.")
            if (length(gMap) < 1 && map == "genetic")
              stop("The genetic map in the input object is empty. Try setting 'map=\"physical\"'")
            if (length(pMap) < 1 && map == "physical")
              stop("The physical map in the input object is empty. Try setting 'map=\"genetic\"'")

            mNames <- mChr <- mLoc <- NULL
            if (length(gMap) > 0 && map == "genetic") {
              mNames <- unlist(sapply(gMap, names), use.names=FALSE)
              mLoc <- unlist(gMap, use.names=FALSE)
              names(mLoc) <- mNames
              elen <- sapply(gMap, length)
              mChr <- rep(names(gMap), times=elen)
              names(mChr) <- mNames
            }

            if (length(pMap) > 0 && map == "physical") {
              mNames <- unlist(sapply(pMap, names), use.names=FALSE)
              mLoc <- unlist(pMap, use.names=FALSE)
              names(mLoc) <- mNames
              elen <- sapply(pMap, length)
              mChr <- rep(names(gMap), times=elen)
              names(mChr) <- mNames
            }

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

setMethod("ciseQTL", signature(x="eQTLnetwork", cisr="numeric"),
          function(x, cisr) {
            eqtl <- alleQTL(x, map="physical")
            eqtl$chrom <- rankSeqlevels(eqtl$chrom)

            ## get the genes and their annotations involved
            eqtlgenes <- geneAnnotation(x)[eqtl$gene]
            mcols(eqtlgenes) <- cbind(mcols(eqtlgenes),
                                      DataFrame(seqnamesRnk=rankSeqlevels(as.vector(seqnames(eqtlgenes)))))

            ## add chromosome and location of the TSS of each eQTL gene
            ciseqtl <- cbind(eqtl,
                             genechrom=eqtlgenes$seqnamesRnk,
                             genelocation=ifelse(strand(eqtlgenes) == "+",
                                                 start(eqtlgenes),
                                                 end(eqtlgenes)))

            ## select those eQTLs occurring on the same chromosome as their target
            ## gene and within a cis radius absolute distance smaller or equal than 'cisr'
            ciseqtl <- ciseqtl[ciseqtl$chrom == ciseqtl$genechrom &
                               abs(ciseqtl$location - ciseqtl$genelocation) <= cisr, ]

            ciseqtl <- ciseqtl[, 1:4]
            rownames(ciseqtl) <- 1:nrow(ciseqtl)
            ciseqtl
          })

setMethod("resetCutoffs", signature(object="eQTLnetwork"),
          function(object) {
            object@p.value <- NA_real_
            object@epsilon <- NA_real_
            object@alpha <- NA_real_
            g.0 <- graphBAM(df=data.frame(from=character(), to=character(), weight=integer()))
            object@qpg <- new("qpGraph", p=0L, q=integer(), n=NA_integer_, epsilon=NA_real_, g=g.0)
            object
          })

## plot method
setMethod("plot", signature(x="eQTLnetwork"),
          function(x, type=c("dot"), xlab="eQTL location", ylab="Gene location", axes=TRUE, ...) {
            type <- match.arg(type)

            if (type == "dot") {
              eqtl <- alleQTL(x, map="physical")

              ## get the genes and their annotations involved
              eqtlgenes <- geneAnnotation(x)[unique(eqtl$gene)]
              mcols(eqtlgenes) <- cbind(mcols(eqtlgenes),
                                        DataFrame(seqnamesRnk=rankSeqlevels(as.vector(seqnames(eqtlgenes)))))

              ## get the loci involved
              qtl <- unique(eqtl[, c("chrom", "location", "QTL")])
              rownames(qtl) <- qtl$QTL
              qtl <- qtl[, c("chrom", "location")]
              qtl$chrom <- rankSeqlevels(qtl$chrom)

              ## re-scale genes and marker locations in eQTL, between 0 and 1
              chrinvolved <- sort(unique(intersect(qtl$chrom, eqtlgenes$seqnamesRnk)))
              chrLen <- seqlengths(eqtlgenes)[chrinvolved]
              chrRelLen <- chrLen / sum(chrLen)
              chrRelCumLen <- c(0, cumsum(chrRelLen)[-length(chrRelLen)])
              geneRelPos <- chrRelCumLen[eqtlgenes$seqnamesRnk] +
                            chrRelLen[eqtlgenes$seqnamesRnk]*(start(eqtlgenes) / chrLen[eqtlgenes$seqnamesRnk])
              names(geneRelPos) <- names(eqtlgenes)

              qtlRelPos <- chrRelCumLen[qtl$chrom] +
                           chrRelLen[qtl$chrom]*(qtl$location/chrLen[qtl$chrom]) 
              names(qtlRelPos) <- rownames(qtl)

              ## perform the dot plot
              plot(qtlRelPos[eqtl$QTL], geneRelPos[eqtl$gene], pch=".",
                   xlim=c(0, 1), ylim=c(0, 1), xlab=xlab, ylab=ylab, cex=4, axes=FALSE, ...)
              segments(c(chrRelCumLen, 1), 0, c(chrRelCumLen, 1), 1, col=gray(0.75), lty="dotted", lwd=2)
              segments(0, c(chrRelCumLen, 1), 1, c(chrRelCumLen, 1), col=gray(0.75), lty="dotted", lwd=2)
              if (axes) {
                axis(1, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
                          labels=names(chrLen), tick=FALSE, cex.axis=0.9)
                axis(2, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
                          labels=names(chrLen), tick=FALSE, cex.axis=0.9, las=1)
              }
            }
          })
