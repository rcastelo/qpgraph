## pretty printing of q orders when there are more than 5, using integer ranges
qPretty <- function(object, cutoff=5) {
  qstr <- ""
  if (!is.na(object@p.value))
    qstr <- "0"

  qstr2 <- paste(object@qOrders, collapse=",")
  if (length(object@qOrders) > 5) {
     dif <- c(object@qOrders[-1], object@qOrders[length(object@qOrders)]+2) - object@qOrders
     qstr2 <- "" ; il <- NA ; ir <- NA
     for (i in seq(along=dif)) {
       if (dif[i] == 1 && is.na(il))
         il <- object@qOrders[i]
       else if (dif[i] > 1 && !is.na(il)) {
         ir <- object@qOrders[i]
         qstr2 <- paste(qstr2, sprintf("%d-%d", il, ir), sep=",")
         il <- ir <- NA
       } else if (dif[i] > 1 && is.na(il))
         qstr2 <- paste(qstr2, object@qOrders[i], sep=",")
     }
     qstr2 <- gsub("^,", "", qstr2)
  }

  if (nchar(qstr) > 0) { ## must be "0"
    qstr <- sprintf("0,%s", qstr2)
  } else {
    qstr <- qstr2
  }

  qstr
}

## eQTLnetwork show method
setMethod("show", signature(object="eQTLnetwork"),
          function(object) {
            cat("eQTLnetwork object:\n")
            if (length(object@organism) > 0)
              cat(sprintf("  Organism: %s\n", object@organism))
            if (length(object@genome) > 0)
              cat(sprintf("  Genome: %s\n", unique(genome(object@genome))))
            if (length(object@geneAnnotationTable) > 0)
              cat(sprintf("  Gene annotation: %s\n", object@geneAnnotationTable))
            mNames <- character()
            if (length(geneticMap(object)) > 0)
              mNames <- unlist(sapply(geneticMap(object), names), use.names=FALSE)
            else if (length(physicalMap(object)) > 0)
              mNames <- unlist(sapply(physicalMap(object), names), use.names=FALSE)
            else
              warning("Both the genetic and the physical map are empty.")

            gNames <- names(object@geneAnnotation)
            cat(sprintf("  Input size: %d markers %d genes\n", length(mNames), length(gNames)))
            fnoq <- gsub("\\(.*\\)", "", paste(deparse(object@modelFormula), collapse=""))
            qstr <- qPretty(object, cutoff=5)
            if (nchar(qstr) > 0)
              cat(sprintf("  Model formula: %s (q = %s)\n", fnoq, qstr))
            else
              cat(sprintf("  Model formula: %s\n", fnoq))

            g <- object@qpg@g
            if (numNodes(g) > 0) {
              edg <- extractFromTo(g)
              edg$from <- as.character(edg$from)
              edg$to <- as.character(edg$to)
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

              connectedVtc <- unique(c(edg$from, edg$to))
              cat(sprintf("\n%sand involving %d genes and %d eQTLs\n",
                          padstr,
                          length(intersect(connectedVtc, gNames)),
                          length(intersect(connectedVtc, mNames))))
            }
          })

## getter and setter methods
setMethod("geneNames", signature(object="eQTLnetwork"),
          function(object) {
            gnames <- character()
            g <- object@qpg@g
            if (numNodes(g) > 0)
              gnames <- intersect(nodes(g), names(geneAnnotation(object)))
            gnames
          })

setMethod("markerNames", signature(object="eQTLnetwork"),
          function(object) {
            mNames <- character()
            if (length(geneticMap(object)) > 0)
              mNames <- unlist(sapply(geneticMap(object), names, simplify=FALSE), use.names=FALSE)
            else if (length(physicalMap(object)) > 0)
              mNames <- unlist(sapply(physicalMap(object), names, simplify=FALSE), use.names=FALSE)
            else
              warning("Both the genetic and the physical map are empty.")

            g <- object@qpg@g
            if (numNodes(g) > 0)
              mNames <- intersect(nodes(g), mNames)

            mNames
          })

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

setMethod("graph", signature(object="eQTLnetwork"),
          function(object) {
            object@qpg@g
          })

setMethod("alleQTL", signature(x="eQTLnetwork"),
          function(x, map=c("genetic", "physical"), gene.loc=FALSE) {
            map <- match.arg(map)
            gMap <- geneticMap(x)
            pMap <- physicalMap(x)

            if (gene.loc && map == "genetic") {
              warning("setting map=\"physical\" because gene.loc=TRUE")
              map <- "physical"
            }

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

            cnames <- c("chrom", "location", "QTL", "gene")
            if (gene.loc)
              cnames <- c(cnames, "genechrom", "genelocation", "distance")

            eqtls <- as.data.frame(matrix(NA, nrow=0, ncol=length(cnames),
                                          dimnames=list(NULL, cnames)))
            n.eqtl <- nrow(mgedg)
            if (n.eqtl > 0) {
              eqtls <- data.frame(chrom=rep(NA_character_, n.eqtl),
                                  location=rep(NA_real_, n.eqtl),
                                  QTL=mgedg$from,
                                  gene=mgedg$to,
                                  stringsAsFactors=FALSE)
              eqtls$chrom <- rankSeqlevels(mChr[eqtls$QTL])
              eqtls$location <- mLoc[eqtls$QTL]
              if (gene.loc) {
                ## get the genes and their annotations involved
                genesGR <- geneAnnotation(x)[eqtls$gene]
                mcols(genesGR) <- cbind(mcols(genesGR),
                                        DataFrame(seqnamesRnk=rankSeqlevels(IRanges::as.vector(seqnames(genesGR)))))
                ## add chromosome and location of the TSS of each eQTL gene
                eqtls <- cbind(eqtls,
                               genechrom=genesGR$seqnamesRnk,
                               genelocation=ifelse(strand(genesGR) == "+",
                                                   start(genesGR),
                                                   end(genesGR)))
                eqtls$distance <- rep(Inf, times=nrow(eqtls))
                masksamechrom <- eqtls$chrom == eqtls$genechrom
                eqtls$distance[masksamechrom] <- abs(eqtls$location[masksamechrom] - eqtls$genelocation[masksamechrom])
              }
            }
            eqtls 
          })

setMethod("ciseQTL", signature(x="eQTLnetwork", cisr="missing"),
          function(x, cisr) {
            ciseQTL(x, cisr=1000L)
          })

setMethod("ciseQTL", signature(x="eQTLnetwork", cisr="numeric"),
          function(x, cisr) {
            eqtls <- alleQTL(x, map="physical", gene.loc=TRUE)
            ciseqtls <- eqtls[eqtls$distance <= cisr, ]
            rownames(ciseqtls) <- 1:nrow(ciseqtls)
            ciseqtls
          })

setMethod("allGeneAssociations", signature(x="eQTLnetwork"),
          function(x, gene.loc=FALSE) {
            gNames <- names(x@geneAnnotation)
            gNodes <- intersect(gNames, nodes(x@qpg@g))
            g <- subGraph(gNodes, x@qpg@g)

            alledg <- extractFromTo(g)
            alledg <- data.frame(genefrom=as.character(alledg$from),
                                 geneto=as.character(alledg$to),
                                 stringsAsFactors=FALSE)

            if (gene.loc) {
              ## get the genes and their annotations involved
              genesGR <- geneAnnotation(x)[gNodes]
              mcols(genesGR) <- cbind(mcols(genesGR),
                                      DataFrame(seqnamesRnk=rankSeqlevels(IRanges::as.vector(seqnames(genesGR)))))
              ## add chromosome and location of the TSS of each gene
              mtfrom <- match(alledg$genefrom, names(genesGR))
              mtto <- match(alledg$geneto, names(genesGR))
              alledg <- cbind(alledg,
                              genefromchrom=genesGR$seqnamesRnk[mtfrom],
                              genefromlocation=ifelse(strand(genesGR)[mtfrom] == "-",
                                                      end(genesGR)[mtfrom],
                                                      start(genesGR)[mtfrom]),
                              genetochrom=genesGR$seqnamesRnk[mtto],
                              genetolocation=ifelse(strand(genesGR)[mtto] == "-",
                                                    end(genesGR)[mtto],
                                                    start(genesGR)[mtto]))
            }

            alledg
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
          function(x, type=c("dot", "hive"), chr=1, xlab="eQTL location", ylab="Gene location", axes=TRUE, map="physical", ...) {
            type <- match.arg(type)

            if (type == "dot") {
              eqtl <- alleQTL(x, map=map)

              ## get the genes and their annotations involved
              eqtlgenes <- geneAnnotation(x)[unique(eqtl$gene)]
              mcols(eqtlgenes) <- cbind(mcols(eqtlgenes),
                                        DataFrame(seqnamesRnk=rankSeqlevels(IRanges::as.vector(seqnames(eqtlgenes)))))

              ## get the loci involved
              qtl <- unique(eqtl[, c("chrom", "location", "QTL")])
              rownames(qtl) <- qtl$QTL
              qtl <- qtl[, c("chrom", "location")]
              qtl$chrom <- rankSeqlevels(as.character(qtl$chrom))

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
            } else if (type == "hive") {
              eqtl <- alleQTL(x, map=map, gene.loc=TRUE)
              eqtl <- eqtl[eqtl$chrom == chr, , drop=FALSE]

              ## split eQTL between cis (same chromosome) and trans (different chromosomes)
              ciseqtl <- eqtl[eqtl$chrom == eqtl$genechrom, c("QTL", "gene")]
              transeqtl <- eqtl[eqtl$chrom != eqtl$genechrom, c("QTL", "gene")]
              gg <- allGeneAssociations(x, gene.loc=TRUE)
              gg <- gg[gg$genefromchrom == chr | gg$genetochrom == chr, , drop=FALSE]
              if (any(gg$genefromchrom != chr))
                gg[gg$genefromchrom != chr, ] <- gg[gg$genefromchrom != chr,
                                                    c("geneto", "genefrom",
                                                      "genetochrom", "genetolocation",
                                                      "genefromchrom", "genefromlocation"), drop=FALSE]
              gg <- gg[, c("genefrom", "geneto")]

              pMap <- physicalMap(x)
              markerPos <- data.frame(chromosome=rep(rankSeqlevels(names(pMap)), times=qtl::nmar(pMap)),
                                      position=unlist(pMap, use.names=FALSE))
              rownames(markerPos) <- unlist(sapply(pMap, names), use.names=FALSE)

              genePos <- data.frame(chromosome=rankSeqlevels(IRanges::as.vector(seqnames(sort(x@geneAnnotation)))),
                                    position=ifelse(strand(sort(x@geneAnnotation)) == "-",
                                                    end(sort(x@geneAnnotation)),
                                                    start(sort(x@geneAnnotation))))

              rownames(genePos) <- names(sort(x@geneAnnotation))

              .ploteQTLnetworkHive(eqtlnet=list(cis=ciseqtl, trans=transeqtl, gg=gg),
                                   mp=markerPos, gp=genePos, cl=seqlengths(x@geneAnnotation),
                                   chr=chr, ...)
            }

          })


## assumes chromosome lengths are correctly ordered by chromosome
.ploteQTLnetworkHive <- function(eqtlnet, mp, gp, cl, chr,
                                 tfIDs=NULL, rbpIDs=NULL, seedGenes=NULL,
                                 axis.cols=c("black", "black", "black"),
                                 def.edges.col=gray(0.8), chr.cols=c("black", "grey"),
                                 size.nodes=0.05, weight.edges=1, ...) {

  cis <- eqtlnet$cis
  trans <- eqtlnet$trans
  gg <- eqtlnet$gg

  ng <- nrow(gp) ## number of annotated genes

  ## calculate gene relative position to the interval [0, 1]
  chrRelLen <- cl/sum(cl)
  chrRelCumLen <- c(0, cumsum(chrRelLen)[-length(chrRelLen)])
  names(chrRelCumLen) <- paste("chr", seq(along=cl), sep="")
  grp <- chrRelCumLen[gp[, 1]] +
         chrRelLen[gp[, 1]] * (gp[, 2]/cl[gp[, 1]])
  names(grp) <- rownames(gp)

  mp1 <- mp[mp[, 1] == chr, ]
  gp1 <- gp[gp[, 1] == chr, ]

  ## add dummy markers at first and last chromosome position
  ## to force the functions below plotting the entire chromosome for
  ## the marker axis
  mp1 <- rbind(data.frame(chromosome=chr, position=0, row.names="dummyBeginning"),
               mp1,
               data.frame(chromosome=chr, position=cl[chr], row.names="dummyEnd"))

  ## add dummy gene annotations at first and last chromosome position
  ## to force the functions below plotting the entire chromosome for
  ## the self-chromosome gene axis
  gp1 <- rbind(data.frame(chromosome=chr, position=0, row.names="dummyBeginning"),
               gp1,
               data.frame(chromosome=chr, position=cl[chr], row.names="dummyEnd"))

  rmp <- mp1[,2]/cl[chr] ## normalize position to the [0-1] scale
  rgp <- gp1[,2]/cl[chr] ## normalize position to the [0-1] scale

  ## calculate relative chromosome positions
  cis1c <- match(cis[ ,1], rownames(mp1))
  cis2c <- match(cis[ ,2], rownames(gp1)) + nrow(mp1)
  trans1c <- match(trans[ ,1], rownames(mp1))
  trans2c <- match(trans[ ,2], names(grp)) + nrow(mp1) + nrow(gp1)

  ## annotate eQTLs by gene functional class (TFs or RBPs) or seed
  ctf <- cis[ ,2] %in% tfIDs
  cistf1 <- cis1c[ctf]
  cistf2 <- cis2c[ctf]
  
  crbp <- cis[ ,2] %in% rbpIDs
  cisrbp1 <- cis1c[crbp]
  cisrbp2 <- cis2c[crbp]

  cseed <- cis[, 2] %in% seedGenes
  cisseed1 <- cis1c[cseed]
  cisseed2 <- cis2c[cseed]
  
  ttf <- trans[ ,2] %in% tfIDs
  transtf1 <- trans1c[ttf]  
  transtf2 <- trans2c[ttf]
  
  trbp <- trans[ ,2] %in% rbpIDs
  transrbp1 <- trans1c[trbp]  
  transrbp2 <- trans2c[trbp]
  
  tseed <- trans[, 2] %in% seedGenes
  transseed1 <- trans1c[tseed]
  transseed2 <- trans2c[tseed]
  
  ## among gene-gene associations select those with cis-eQTLs in this
  ## chromosome and annotate function or seed
  gene1_cis <- gene2_cis <- gene1_ctf <- gene2_ctf <- gene1_crbp <- gene2_crbp <- gene1_seed <- gene2_seed <- c()
  cisg <- gg[, 1] %in% cis[, 2]
  cisg2 <- gg[, 2] %in% cis[, 2]
  genesCis <- matrix(NA, nrow=0, ncol=ncol(gg))
  if (sum(cisg) > 0) {
    genesCis <- rbind(gg[cisg, , drop=FALSE], gg[cisg2, 2:1, drop=FALSE])
    gene1_cis <- match(genesCis[, 1], rownames(gp1)) + nrow(mp1)
    gene2_cis <- match(genesCis[, 2], names(grp)) + nrow(mp1) + nrow(gp1)
      
    idtfg <- (genesCis[, 1] %in% tfIDs & !(genesCis[, 1] %in% seedGenes)) |
             (genesCis[, 2] %in% tfIDs & !(genesCis[, 2] %in% seedGenes))
    tfg <- genesCis[idtfg, ,drop=FALSE]      
    gene1_ctf <- match(tfg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_ctf <- match(tfg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)
  
    idrbpg <- (genesCis[, 1] %in% rbpIDs & !(genesCis[, 1] %in% seedGenes)) |
              (genesCis[, 2] %in% rbpIDs & !(genesCis[, 2] %in% seedGenes))
    rbpg <- genesCis[idrbpg, ,drop=FALSE]      
    gene1_crbp <- match(rbpg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_crbp <- match(rbpg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)  

    idseedg <- (genesCis[, 1] %in% seedGenes | genesCis[, 2] %in% seedGenes)
    seedg <- genesCis[idseedg, ,drop=FALSE]      
    gene1_cseed <- match(seedg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_cseed <- match(seedg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)  
  }

  ## among gene-gene associations select those with trans-eQTLs in other
  ## chromosomes and annotate function
  gene1_trans <- gene2_trans <- gene1_ttf <- gene2_ttf <- gene1_trbp <- gene2_trbp <- gene1_seed <- gene2_seed <- gene1_tseed <- gene2_tseed <- c()
  transg <- gg[, 2] %in% trans[, 2]
  transg2 <- gg[, 1] %in% trans[, 2]
  genesTrans <- matrix(NA, nrow=0, ncol=ncol(gg))
  if (sum(transg) > 0) {
    genesTrans <- rbind(gg[transg, , drop=FALSE], gg[transg2, 2:1, drop=FALSE])
    gene1_trans <- match(genesTrans[, 1], rownames(gp1)) + nrow(mp1)
    gene2_trans <- match(genesTrans[, 2], names(grp)) + nrow(mp1) + nrow(gp1)
    
    idtfg <- (genesTrans[, 1] %in% tfIDs & !(genesTrans[, 1] %in% seedGenes)) |
             (genesTrans[, 2] %in% tfIDs & !(genesTrans[, 2] %in% seedGenes))
    tfg <- genesTrans[idtfg, ,drop=FALSE]      
    gene1_ttf <- match(tfg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_ttf <- match(tfg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)

    idrbpg <- (genesTrans[, 1] %in% rbpIDs & !(genesTrans[, 1] %in% seedGenes)) |
              (genesTrans[, 2] %in% rbpIDs & !(genesTrans[, 2] %in% seedGenes))
    rbpg <- genesTrans[idrbpg, ,drop=FALSE]      
    gene1_trbp <- match(rbpg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_trbp <- match(rbpg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)  

    idseedg <- (genesTrans[, 1] %in% seedGenes | genesTrans[, 2] %in% seedGenes)
    seedg <- genesTrans[idseedg, ,drop=FALSE]      
    gene1_tseed <- match(seedg[, 1], rownames(gp1)) + nrow(mp1)
    gene2_tseed <- match(seedg[, 2], names(grp)) + nrow(mp1) + nrow(gp1)  
  }
  
  ## create a list object with all the information to be plotted with
  ## the function '.eQTLplotHive()' that has been adapted from the
  ## R/CRAN package HiveR
  col_nodes <- rep(chr.cols[1], nrow(mp1) + nrow(gp1) + ng) 
  col_nodes[c(rep(FALSE, nrow(mp1)+ nrow(gp1)), (gp[ ,1] %% 2) == 0)] <- chr.cols[2]
  size_nodes <- rep(size.nodes, nrow(mp1) + nrow(gp1) + ng)
    
  id1 <- c(cis1c, trans1c, gene1_cis, gene1_trans,
           cistf1, transtf1, gene1_ctf, gene1_ttf,
           cisrbp1, transrbp1, gene1_crbp, gene1_trbp,
           cisseed1, transseed1, gene1_cseed, gene1_tseed)
  id2 <- c(cis2c, trans2c, gene2_cis, gene2_trans,
           cistf2, transtf2, gene2_ctf, gene2_ttf,
           cisrbp2, transrbp2, gene2_crbp, gene2_trbp,
           cisseed2, transseed2, gene2_cseed, gene2_tseed)
  
  col_edges <- c(rep(def.edges.col, length(c(cis1c, trans1c, gene1_cis, gene1_trans))), 
                 rep("orange", length(c(cistf1, transtf1, gene1_ctf, gene1_ttf))),  
                 rep("blue", length(c(cisrbp1, transrbp1, gene1_crbp, gene1_trbp))),
                 rep("black", length(c(cisseed1, transseed1, gene1_cseed, gene1_tseed))))   


  nodes <- data.frame(id=1:(nrow(mp1) + nrow(gp1) + ng), 
                      lab=c(rownames(mp1), rownames(gp1), names(grp)),
                      axis=as.integer(c(rep(1, nrow(mp1)), rep(2, nrow(gp1)), rep(3, ng))), 
                      radius=c(rmp, rgp, grp),
                      size=size_nodes, 
                      color=col_nodes,
                      stringsAsFactors=FALSE)

  edges <- data.frame(id1=id1,
      	              id2=id2,
                      weight=rep(weight.edges, length(id1)), 
                      color=col_edges, 
                      stringsAsFactors=FALSE)

  hpd <- list(nodes=nodes,       
              edges=edges,
              type="2D",
              desc="eQTL Hive Plot",
              axis.cols=axis.cols)
  class(hpd) <- "HivePlotData"
  ## stopifnot(!chkHPD(hpd))

  .eQTLplotHive(hpd, ch=0.1, dr.nodes=TRUE, np=FALSE, ...)
}

## adapted from the plotHive function in the R/CRAN package HiveR by Bryan Hanson
## see http://cran.r-project.org/web/packages/HiveR/index.html
.eQTLplotHive <- function(HPD, ch = 1, ## method = "abs", ## MODIFICATION WE DON'T NEED THIS ARGUMENT
	dr.nodes = TRUE, bkgnd = "white",
	axLabs = NULL, axLab.pos = NULL, axLab.gpar = NULL,
	anNodes = NULL, anNode.gpar = NULL,
	arrow = NULL, np = TRUE, ...) {
	
	# Function to plot hive plots using grid graphics
	# Inspired by the work of Martin Kryzwinski
	# Bryan Hanson, DePauw Univ, Feb 2011 onward
	
	# This function is intended to draw in 2D for nx from 2 to 6
	# The results will be similar to the original hive plot concept

##### Set up some common parameters

	if (!HPD$type == "2D") stop("This is not a 2D hive data set: use plot3dHive instead from package HiveR")
	## chkHPD(HPD)
	nx <- length(unique(HPD$nodes$axis))

	if (nx == 1) stop("Something is wrong: only one axis seems to be present")

	# Send out for ranking/norming/pruning/inverting if requested
	
  ## MODIFICATION WE DON'T NEED THIS ARGUMENT, IN PRINCIPLE
	## if (!method == "abs") HPD <- manipAxis(HPD, method, ...)

	nodes <- HPD$nodes
	edges <- HPD$edges
	axis.cols <- HPD$axis.cols

	# Fix up center hole
	
	nodes$radius <- nodes$radius + ch
	HPD$nodes$radius <- nodes$radius

##### Some convenience functions, only defined in this function environ.
##### The two long functions need to stay here for simplicity, since
##### all of the radius checking etc is here and if moved elsewhere,
##### these calculations would have to be redone or results passed.

	p2cX <- function(r, theta) x <- r*cos(theta*2*pi/360)
	p2cY <- function(r, theta) y <- r*sin(theta*2*pi/360)

	addArrow <- function(arrow, nx) {
		if (!length(arrow) >= 5) stop("Too few arrow components")
		if (is.null(axLab.gpar)) {
			if (bkgnd == "black") axLab.gpar <- gpar(fontsize = 12, col = "white", lwd = 2)
			if (!bkgnd == "black") axLab.gpar <- gpar(fontsize = 12, col = "black", lwd = 2)
			}
		a <- as.numeric(arrow[2])
		rs <- as.numeric(arrow[3])
		re <- as.numeric(arrow[4])
		b <- as.numeric(arrow[5]) # label offset from end of arrow
		
		x.st <- p2cX(rs, a)
		y.st <- p2cY(rs, a)
		x.end <- p2cX(re, a)
		y.end <- p2cY(re, a)
					
		x.lab <- p2cX(re + b, a) # figure arrow label position
		y.lab <- p2cY(re + b, a)
		al <- 0.2*(re-rs) # arrow head length
		
		# for nx = 2 only, offset the arrow
		# in the y direction to save space overall
		
		if (nx == 2) {
			if (is.na(arrow[6])) {
				arrow[6] <- 0
				cat("\tThe arrow can be offset vertically; see ?plotHive\n")
				}
			y.st <- y.st + as.numeric(arrow[6])
			y.end <- y.end + as.numeric(arrow[6])
			y.lab <- y.lab + as.numeric(arrow[6])		
			}

		grid.lines(x = c(x.st, x.end), y = c(y.st, y.end),
			arrow = arrow(length = unit(al, "native")),
			default.units = "native", gp = axLab.gpar)
		grid.text(arrow[1], x.lab, y.lab, default.units = "native", gp = axLab.gpar)
		}

	annotateNodes <- function(anNodes, nodes, nx) {

		if (is.null(anNode.gpar)) {
			if (bkgnd == "black") anNode.gpar <- gpar(fontsize = 10, col = "white", lwd = 0.5)
			if (!bkgnd == "black") anNode.gpar <- gpar(fontsize = 10, col = "black", lwd = 0.5)
			}
			
    ## MODIFICATION for eQTL plotting purposes
    ## do not read gene annotations from CSV files but
    ## directly from a data.frame object
		## ann <- read.csv(anNodes, header = TRUE)
		# Columns should be: node.lab, node.text, angle, radius, offset
		# Locate the node on an axis at a particular radius
    ann <- anNodes
		
		id <- c()	
		for (n in 1:nrow(ann)) {			
			pat <- paste("\\b", ann$node.lab[n], "\\b", sep = "") 
			## id <- c(id, grep(pat, nodes$lab))
			id <- c(id, grep(pat, nodes$lab)[1]) ## ROBERT: take only the first axis labels
			}
		
			N <- matrix(data = c(
			0, 180, NA, NA, NA, NA,
			90, 210, 330, NA, NA, NA,
			90, 180, 270, 0, NA, NA,
			90, 162, 234, 306, 18, NA,
			90, 150, 210, 270, 330, 390),
			byrow = TRUE, nrow = 5)
			
		ax <- nodes$axis[id] # axis number
		for (n in 1:length(ax)) {
			ax[n] <- N[nx-1,ax[n]]			
			}

		x.st <- p2cX(nodes$radius[id], ax)
		y.st <- p2cY(nodes$radius[id], ax)
		
		x.end <- p2cX(ann$radius, ann$angle)
		y.end <- p2cY(ann$radius, ann$angle)
					
		x.lab <- p2cX(ann$radius + ann$offset, ann$angle) # figure label position
		y.lab <- p2cY(ann$radius + ann$offset, ann$angle)
		
    ## MODIFICATION for eQTL plotting purposes
    ## do not annotate with segments, draw text right on the axis
		## grid.segments(x0 = x.st, x1 = x.end, y0 = y.st, y1 = y.end,
		## default.units = "native", gp = anNode.gpar)
		## grid.text(ann$node.text, x.lab, y.lab, hjust = ann$hjust, vjust = ann$vjust,
		grid.text(ann$node.text, x.st, y.st, hjust = ann$hjust, vjust = ann$vjust,
			default.units = "native", gp = anNode.gpar)
		}

###############

	# Figure out which nodes to draw for each edge
	# Since they are in random order
	# Do this once/early to save time
	
	id1 <- id2 <- c()
	
	for (n in 1:nrow(edges)) {
		
		pat1 <- paste("\\b", edges$id1[n], "\\b", sep = "") 
		pat2 <- paste("\\b", edges$id2[n], "\\b", sep = "")
		id1 <- c(id1, grep(pat1, nodes$id))
		id2 <- c(id2, grep(pat2, nodes$id))
		
		}

##### Two dimensional case  (using grid graphics)

	# Prep axes first
	
	if (nx == 2) {
		
		n1 <- subset(nodes, axis == 1)
		n2 <- subset(nodes, axis == 2)
		max1 <- max(n1$radius)
		max2 <- max(n2$radius)
		min1 <- min(n1$radius)
		min2 <- min(n2$radius)
	
		r.st <- c(min1, min2) # in polar coordinates
		axst <- c(0, 180)
		x0a = p2cX(r.st, axst)
		y0a = p2cY(r.st, axst)

		r.end <- c(max1, max2)
		axend <- c(0, 180)
		x1a = p2cX(r.end, axend)
		y1a = p2cY(r.end, axend)
	
	# Set up grid graphics viewport
		
		md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.5 # max dimension
		# 1.5 is used in case of labels
		
		if (np) grid.newpage()
		grid.rect(gp = gpar(fill = bkgnd))
		vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
			xscale = c(-md, md), yscale = c(-md, md),
			name = "3DHivePlot")

		pushViewport(vp)
	
	# Now draw edges
		
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if (nodes$axis[id1[n]] == 1) { # set up edge start params 1st
				th.st <- c(th.st, 0)
				r.st <- c(r.st, nodes$radius[id1[n]])
				}
			if (nodes$axis[id1[n]] == 2) {
				th.st <- c(th.st, 180)
				r.st <- c(r.st, nodes$radius[id1[n]])
				}

			if (nodes$axis[id2[n]] == 1) { # now edge end params
				th.end <- c(th.end, 0)
				r.end <- c(r.end, nodes$radius[id2[n]])
				}
			if (nodes$axis[id2[n]] == 2) {
				th.end <- c(th.end, 180)
				r.end <- c(r.end, nodes$radius[id2[n]])
				}

			ecol <- c(ecol, edges$color[n])
			ewt <- c(ewt, edges$weight[n])
			}
				
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Draw axes
	
		grid.segments(x0a, y0a, x1a, y1a,
			gp = gpar(col = HPD$axis.cols, lwd = 8),
			default.units = "native")
	
	# Now add nodes
	
		if (dr.nodes) {
			r <- c(n1$radius, n2$radius) 
			theta <- c(rep(0, length(n1$radius)),
				rep(180, length(n2$radius)))
			x = p2cX(r, theta)
			y = p2cY(r, theta)
			grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size), col = c(n1$color, n2$color)))
			}

	# Now label axes
		
		if (!is.null(axLabs)) {
			if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
			if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
			r <- c(max1, max2)
			if (is.null(axLab.pos)) axLab.pos <- r*0.1
			r <- r + axLab.pos
			th <- c(0, 180)
			x <- p2cX(r, th)
			y <- p2cY(r, th)
			grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
			}

	# Add a legend arrow & any annotations
		
		if (!is.null(arrow)) addArrow(arrow, nx)
		if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx)

		} # end of 2D
	
##### Three dimensional case (using grid graphics)

	# Prep axes first
	
	if (nx == 3) {
		
		n1 <- subset(nodes, axis == 1)
		n2 <- subset(nodes, axis == 2)
		n3 <- subset(nodes, axis == 3)
		max1 <- max(n1$radius)
		max2 <- max(n2$radius)
		max3 <- max(n3$radius)
		min1 <- min(n1$radius)
		min2 <- min(n2$radius)
		min3 <- min(n3$radius)

		r.st <- c(min1, min2, min3) # in polar coordinates
		axst <- c(90, 210, 330)
		x0a = p2cX(r.st, axst)
		y0a = p2cY(r.st, axst)

		r.end <- c(max1, max2, max3)
		axend <- c(90, 210, 330)
		x1a = p2cX(r.end, axend)
		y1a = p2cY(r.end, axend)
	
	# Set up grid graphics viewport
	
		md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
		if (np) grid.newpage()
		grid.rect(gp = gpar(fill = bkgnd))
		vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
			xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
		pushViewport(vp)

	# Now draw edges (must do in sets as curvature is not vectorized)
		
	# Axis 1 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 2 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 3 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 1 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 3 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 2 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Draw axes
	
		grid.segments(x0a, y0a, x1a, y1a,
			gp = gpar(col = HPD$axis.cols, lwd = 3),
			default.units = "native")
	
	# Now add nodes
	
		if (dr.nodes) {
			r <- c(n1$radius, n2$radius, n3$radius) 
			theta <- c(rep(90, length(n1$radius)),
				rep(210, length(n2$radius)),
				rep(330, length(n3$radius)))
			x = p2cX(r, theta)
			y = p2cY(r, theta)
			grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size),
			col = c(n1$color, n2$color, n3$color)))
			}

	# Now label axes
		
		if (!is.null(axLabs)) {
			if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
			if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
			r <- c(max1, max2, max3)
			if (is.null(axLab.pos)) axLab.pos <- r*0.1
			r <- r + axLab.pos
			th <- c(90, 210, 330)
			x <- p2cX(r, th)
			y <- p2cY(r, th)
			grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
			}

	# Add a legend arrow & any annotations
		
		if (!is.null(arrow)) addArrow(arrow, nx)
		if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx)

		} # end of 3D
	

##### Four dimensional case (using grid graphics)

	# Prep axes first
	
	if (nx == 4) {
		
		n1 <- subset(nodes, axis == 1)
		n2 <- subset(nodes, axis == 2)
		n3 <- subset(nodes, axis == 3)
		n4 <- subset(nodes, axis == 4)
		max1 <- max(n1$radius)
		max2 <- max(n2$radius)
		max3 <- max(n3$radius)
		max4 <- max(n4$radius)
		min1 <- min(n1$radius)
		min2 <- min(n2$radius)
		min3 <- min(n3$radius)
		min4 <- min(n4$radius)

		r.st <- c(min1, min2, min3, min4) # in polar coordinates
		axst <- c(90, 180, 270, 0)
		x0a = p2cX(r.st, axst)
		y0a = p2cY(r.st, axst)

		r.end <- c(max1, max2, max3, max4)
		axend <- c(90, 180, 270, 0)
		x1a = p2cX(r.end, axend)
		y1a = p2cY(r.end, axend)
	
	# Set up grid graphics viewport
	
		md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.5 # max dimension
		if (np) grid.newpage()
		grid.rect(gp = gpar(fill = bkgnd))
		vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
			xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
		pushViewport(vp)

	# Now draw edges (must do in sets as curvature is not vectorized)
		
	# Axis 1 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 180)
				r.end <- c(r.end, nodes$radius[id2[n]])
			ecol <- c(ecol, edges$color[n])
			ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 2 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 180)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
			ecol <- c(ecol, edges$color[n])
			ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 3 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 0)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 4 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 0)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
			ecol <- c(ecol, edges$color[n])
			ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 1 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 0)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 4 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 0)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 3 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 180)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}				
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 2 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 180)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}				
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 180)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 180)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 0)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 0)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Draw axes
	
		grid.segments(x0a, y0a, x1a, y1a,
			gp = gpar(col = HPD$axis.cols, lwd = 3),
			default.units = "native")
	
	# Now add nodes
	
		if (dr.nodes) {
			r <- c(n1$radius, n2$radius, n3$radius, n4$radius) 
			theta <- c(rep(90, length(n1$radius)),
				rep(180, length(n2$radius)),
				rep(270, length(n3$radius)),
				rep(0, length(n4$radius)))
			x = p2cX(r, theta)
			y = p2cY(r, theta)
			grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size),
			col = c(n1$color, n2$color, n3$color, n4$color)))
			}

	# Now label axes
		
		if (!is.null(axLabs)) {
			if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
			if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
			r <- c(max1, max2, max3, max4)
			if (is.null(axLab.pos)) axLab.pos <- r*0.1
			r <- r + axLab.pos
			th <- c(90, 180, 270, 0)
			x <- p2cX(r, th)
			y <- p2cY(r, th)
			grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
			}

	# Add a legend arrow & any annotations
		
		if (!is.null(arrow)) addArrow(arrow, nx)
		if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx)

		} # end of 4D
	
##### Five dimensional case (using grid graphics)

	# Prep axes first
	
	if (nx == 5) {
		
		n1 <- subset(nodes, axis == 1)
		n2 <- subset(nodes, axis == 2)
		n3 <- subset(nodes, axis == 3)
		n4 <- subset(nodes, axis == 4)
		n5 <- subset(nodes, axis == 5)
		max1 <- max(n1$radius)
		max2 <- max(n2$radius)
		max3 <- max(n3$radius)
		max4 <- max(n4$radius)
		max5 <- max(n5$radius)
		min1 <- min(n1$radius)
		min2 <- min(n2$radius)
		min3 <- min(n3$radius)
		min4 <- min(n4$radius)
		min5 <- min(n5$radius)

		r.st <- c(min1, min2, min3, min4, min5) # in polar coordinates
		axst <- c(90, 162, 234, 306, 18)
		x0a = p2cX(r.st, axst)
		y0a = p2cY(r.st, axst)

		r.end <- c(max1, max2, max3, max4, max5)
		axend <- c(90, 162, 234, 306, 18)
		x1a = p2cX(r.end, axend)
		y1a = p2cY(r.end, axend)
	
	# Set up grid graphics viewport
	
		md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
		if (np) grid.newpage()
		grid.rect(gp = gpar(fill = bkgnd))
		vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
			xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
		pushViewport(vp)

	# Now draw edges (must do in sets as curvature is not vectorized)
		
	# Axis 1 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 162)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}			
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 2 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 162)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 234)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 3 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 234)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 306)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 4 -> 5
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 306)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 18)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 5 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 18)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}


	# Axis 1 -> 5
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 18)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 5 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 18)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 306)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 4 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 306)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 234)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 3 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 234)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 162)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 2 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 162)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 162)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 162)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 234)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 234)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 306)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 306)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 18)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 18)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Draw axes
	
		grid.segments(x0a, y0a, x1a, y1a,
			gp = gpar(col = HPD$axis.cols, lwd = 3),
			default.units = "native")
	
	# Now add nodes
	
		if (dr.nodes) {
			r <- c(n1$radius, n2$radius, n3$radius, n4$radius, n5$radius) 
			theta <- c(rep(90, length(n1$radius)),
				rep(162, length(n2$radius)),
				rep(234, length(n3$radius)),
				rep(306, length(n4$radius)),
				rep(18, length(n5$radius)))
			x = p2cX(r, theta)
			y = p2cY(r, theta)
			grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size, n5$size),
			col = c(n1$color, n2$color, n3$color, n4$color, n5$color)))
			}

	# Now label axes
		
		if (!is.null(axLabs)) {
			if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
			if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
			r <- c(max1, max2, max3, max4, max5)
			if (is.null(axLab.pos)) axLab.pos <- r*0.1
			r <- r + axLab.pos
			th <- c(90, 162, 234, 306, 18)
			x <- p2cX(r, th)
			y <- p2cY(r, th)
			grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
			}

	# Add a legend arrow & any annotations
		
		if (!is.null(arrow)) addArrow(arrow, nx)
		if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx)

		} # end of 5D

##### Six dimensional case (using grid graphics)

	# Prep axes first
	
	if (nx == 6) {
		
		n1 <- subset(nodes, axis == 1)
		n2 <- subset(nodes, axis == 2)
		n3 <- subset(nodes, axis == 3)
		n4 <- subset(nodes, axis == 4)
		n5 <- subset(nodes, axis == 5)
		n6 <- subset(nodes, axis == 6)
		max1 <- max(n1$radius)
		max2 <- max(n2$radius)
		max3 <- max(n3$radius)
		max4 <- max(n4$radius)
		max5 <- max(n5$radius)
		max6 <- max(n6$radius)
		min1 <- min(n1$radius)
		min2 <- min(n2$radius)
		min3 <- min(n3$radius)
		min4 <- min(n4$radius)
		min5 <- min(n5$radius)
		min6 <- min(n6$radius)

		r.st <- c(min1, min2, min3, min4, min5, min6) # in polar coordinates
		axst <- c(90, 150, 210, 270, 330, 390)
		x0a = p2cX(r.st, axst)
		y0a = p2cY(r.st, axst)

		r.end <- c(max1, max2, max3, max4, max5, max6)
		axend <- c(90, 150, 210, 270, 330, 390)
		x1a = p2cX(r.end, axend)
		y1a = p2cY(r.end, axend)
	
	# Set up grid graphics viewport
	
		md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
		if (np) grid.newpage()
		grid.rect(gp = gpar(fill = bkgnd))
		vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
			xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
		pushViewport(vp)

	# Now draw edges (must do in sets as curvature is not vectorized)
		
	# Axis 1 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {

			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 150)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 2 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 150)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 3 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 4 -> 5
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 5 -> 6
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 6)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 390)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 6 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 390)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Axis 1 -> 6
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 6)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 390)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 6 -> 5
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 390)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 5 -> 4
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 4 -> 3
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)
		
		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 3 -> 2
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 150)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 2 -> 1
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
						
			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 150)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])				}
			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
			}

	# Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
	
		r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
			
		for (n in 1:nrow(edges)) {
			
			if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
				th.st <- c(th.st, 90)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 90)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
				th.st <- c(th.st, 150)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 150)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
				th.st <- c(th.st, 210)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 210)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}
				
			if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
				th.st <- c(th.st, 270)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 270)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 5)) {
				th.st <- c(th.st, 330)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 330)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 6)) {
				th.st <- c(th.st, 390)
				r.st <- c(r.st, nodes$radius[id1[n]])
				th.end <- c(th.end, 390)
				r.end <- c(r.end, nodes$radius[id2[n]])
				ecol <- c(ecol, edges$color[n])
				ewt <- c(ewt, edges$weight[n])
				}

			}
		
		x0 = p2cX(r.st, th.st)
		y0 = p2cY(r.st, th.st)
		x1 = p2cX(r.end, th.end)
		y1 = p2cY(r.end, th.end)

		if (!length(x0) == 0) {
			grid.curve(x0, y0, x1, y1,
				default.units = "native", ncp = 5, square = FALSE,
				gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
			}

	# Draw axes
	
		# grid.segments(x0a, y0a, x1a, y1a,
			# gp = gpar(col = "black", lwd = 7),
			# default.units = "native") # more like linnet

		grid.segments(x0a, y0a, x1a, y1a,
			gp = gpar(col = HPD$axis.cols, lwd = 3),
			default.units = "native")
	
	# Now add nodes
	
		if (dr.nodes) {
			r <- c(n1$radius, n2$radius, n3$radius, n4$radius, n5$radius, n6$radius) 
			theta <- c(rep(90, length(n1$radius)),
				rep(150, length(n2$radius)),
				rep(210, length(n3$radius)),
				rep(270, length(n4$radius)),
				rep(330, length(n5$radius)),
				rep(390, length(n6$radius)))
			x = p2cX(r, theta)
			y = p2cY(r, theta)
			grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size, n5$size, n6$size),
			col = c(n1$color, n2$color, n3$color, n4$color, n5$color, n6$color)))
			}

	# Now label axes
		
		if (!is.null(axLabs)) {
			if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
			if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
			r <- c(max1, max2, max3, max4, max5, max6)
			if (is.null(axLab.pos)) axLab.pos <- r*0.1
			r <- r + axLab.pos
			th <- c(90, 150, 210, 270, 330, 390)
			x <- p2cX(r, th)
			y <- p2cY(r, th)
			grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
			}

	# Add a legend arrow & any annotations
		
		if (!is.null(arrow)) addArrow(arrow, nx)
		if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx)

		} # end of 6D

	
	} # closing brace, this is the end!
