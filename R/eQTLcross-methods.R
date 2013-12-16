setValidity("eQTLcross",
              ## order vertex alphabetically within each edge (needed to avoid graph-specific nightmare)
            function(object) {
              valid <- TRUE

              if (object@type == "bc" && object@model@dLevels != 2)
                valid <- "when 'type'=\"bc\" the 'HMgmm' object should have 2 levels for the discrete variables."

              valid
            })

## constructor methods: 'map' is a genetic map as defined in the 'qtl' package
setMethod("eQTLcross", signature(map="map", genes="missing", model="missing"),
          function(map, genes, model, type="bc",
                   geneNetwork=matrix(NA, nrow=0, ncol=2),
                   rho=0.5, sigma=diag(0)) {
            type <- match.arg(type)
            genes <- matrix(NA, nrow=0, ncol=2)
            model <- matrix(NA, nrow=0, ncol=4)
            eQTLcross(map, genes, model, type, geneNetwork, rho, sigma)
          })

## model is a matrix where each row corresponds to a different eQTL association and
## contains the chromosome number, cM position, gene and additive effect of the eQTL
## in this constructor both 'model' and 'geneNetwork' should not contain gene names
## but indexes to the rows of the 'genes' matrix describing the genetic map of genes
## whenever the 'genes' matrix contains rownames, those will be used as gene names,
## otherwise they'll be automatically generated. For QTLs their names will be always
## automatically generated as QTL1, QTL2, etc. following their order along the genome.
## this may change in the future to accept arbitrary names for QTLs.
## this constructor stores the structure of eQTL and gene-gene associations in a graphBAM object
setMethod("eQTLcross", signature(map="map", genes="matrix", model="matrix"),
          function(map, genes, model, type="bc",
                   geneNetwork=matrix(NA, nrow=0, ncol=2),
                   rho=0.5, sigma=diag(nrow(genes))) {
            type <- match.arg(type)

            nGenes <- nrow(genes)

            dVertexLabels <- c()
            cVertexLabels <- rownames(genes)
            if (nGenes > 0) {
              if (is.null(rownames(genes))) {
                cVertexLabels <- paste0("g", seq_len(nGenes))
                rownames(genes) <- cVertexLabels
              }
              rownames(sigma) <- colnames(sigma) <- cVertexLabels
            }

            edges <- as.data.frame(matrix(NA, nrow=0, ncol=3,
                                          dimnames=list(NULL, c("from", "to", "weight"))),
                                   stringsAsFactors=FALSE)
            if (nrow(model) > 0) {
              model <- model[order(model[, 1], model[, 2]), , drop=FALSE]             ## sort eQTL by chromosome and location
              rownames(model) <- NULL                                                 ## remove rownames, if any
              eQTLtag <- rle(apply(model[, 1:2, drop=FALSE], 1, paste, collapse="_")) ## tag every different eQTL location
              nQTL <- length(eQTLtag$values)                                          ## count the number of different QTL
              model <- cbind(model, rep(1:nQTL, times=eQTLtag$lengths))               ## add QTL number to each different QTL
              dVertexLabels <- paste0("QTL", 1:nQTL)                                  ## assign QTL names to each different QTL
              edges <- data.frame(from=paste0("QTL", model[, 5]),
                                  to=rownames(genes)[model[, 3]],
                                  weight=1, stringsAsFactors=FALSE)
            }
            
            if (nrow(geneNetwork) > 0)
              edges <- rbind(edges,
                             data.frame(from=rownames(genes)[geneNetwork[, 1]],
                                        to=rownames(genes)[geneNetwork[, 2]],
                                        weight=rep(1, nrow(geneNetwork)),
                                        stringsAsFactors=FALSE))

            g <- graph::graphBAM(edges, nodes=c(dVertexLabels, cVertexLabels))
            if (length(c(dVertexLabels, cVertexLabels)) > 0) {
              graph::nodeDataDefaults(g, "type") <- "continuous"
              graph::nodeDataDefaults(g, "chr") <- NA_integer_
              graph::nodeDataDefaults(g, "loc") <- NA_real_
              graph::edgeDataDefaults(g, "a") <- NA_real_
            }

            if (nrow(model) > 0) {
              qtlchr <- unique(model[, c(1, 5), drop=FALSE])
              qtlloc <- unique(model[, c(2, 5), drop=FALSE])
              qtlchr <- do.call("names<-", list(qtlchr[, 1], paste0("QTL", qtlchr[, 2])))
              qtlloc <- do.call("names<-", list(qtlloc[, 1], paste0("QTL", qtlloc[, 2])))
              graph::nodeData(g, dVertexLabels, "type") <- "discrete"
              graph::nodeData(g, dVertexLabels, "chr") <- qtlchr[dVertexLabels]
              graph::nodeData(g, dVertexLabels, "loc") <- qtlloc[dVertexLabels]
              if (nrow(model) > 0) {
                mixedEdgeMask <- !is.na(match(edges$from, dVertexLabels))
                graph::edgeData(g, from=edges[mixedEdgeMask, ]$from, to=edges[mixedEdgeMask, ]$to, "a") <- model[, 4]
              }
            }

            eQTLcross(map=map, genes=genes, model=g, type=type, rho=rho, sigma=sigma)
          })

setMethod("eQTLcross", signature(map="map", genes="matrix", model="graphBAM"),
          function(map, genes, model,
                   type="bc",
                   rho=0.5,
                   sigma=diag(ifelse(nrow(genes) > 0, sum(unlist(graph::nodeData(model, graph::nodes(model), attr="type"), use.names=FALSE) == "continuous"), 0))) {
            type <- match.arg(type)

            dLevels <- switch(type,
                              bc=2L)
            a <- numeric()
            nGenes <- nrow(genes)

            if (nGenes > 0) {
              maskcvertex <- unlist(graph::nodeData(model, graph::nodes(model), attr="type"), use.names=FALSE) == "continuous"
              cVertexLabels <- graph::nodes(model)[maskcvertex]
              maskpositivedeg <- graph::degree(model, cVertexLabels) > 0

              ## additive effects per gene are the sum of additive effects from each eQTL
              a <- do.call("names<-", list(rep(0, length(cVertexLabels)), cVertexLabels))
              a[maskpositivedeg] <- sapply(cVertexLabels[maskpositivedeg],
                                           function(v) sum(unlist(edgeData(model, to=v, attr="a"), use.names=FALSE),
                                                           na.rm=TRUE))

              if (is.null(rownames(sigma)) || is.null(colnames(sigma)))
                rownames(sigma) <- colnames(sigma) <- cVertexLabels
            }

            gmm <- HMgmm(model, dLevels=dLevels, a=a, rho=rho, sigma=sigma)
            eQTLcross(map, genes=genes, model=gmm, type=type)
          })

setMethod("eQTLcross", signature(map="map", genes="missing", model="HMgmm"),
          function(map, genes, model, type) {

            eQTLcross(map, matrix(NA, nrow=0, ncol=2), model, type)
          })

setMethod("eQTLcross", signature(map="map", genes="matrix", model="HMgmm"),
          function(map, genes, model, type) {

            new("eQTLcross", map=map, genes=genes, model=model, type=type)
          })

setMethod("addGenes", signature(eQTLcross="eQTLcross", n="missing"),
          function(eQTLcross, n, chr=NULL, location=NULL) {
            addGenes(eQTLcross, 1L, chr, location)
          })

setMethod("addGenes", signature(eQTLcross="eQTLcross", n="numeric"),
          function(eQTLcross, n=1, chr=NULL, location=NULL) {
            addGenes(eQTLcross, as.integer(n), chr, location)
          })

setMethod("addGenes", signature(eQTLcross="eQTLcross", n="integer"),
          function(eQTLcross, n=1L, chr=NULL, location=NULL) {
            if (length(n) > 1)
              stop("Second argument 'n' should contain one single value specifying the number of genes to add.")

            if (is.null(chr))
              chr <- sample(seq(along=qtl::chrlen(eQTLcross@map)), size=n, replace=TRUE)
            else if (length(chr) == 1)
              chr <- rep(chr, n)
            else if (length(chr) != n)
              stop("Chromosomes in 'chr' should match the number of genes 'n'.")

            if (is.null(location))
              location <- sapply(chr, function(chr) runif(1, min=0, max=qtl::chrlen(eQTLcross@map)[chr]))
            else if (length(location) != n)
              stop("Chromosomal locations in 'location' should match the number of genes 'n'.")

            if (any(duplicated(location)))
              stop("Duplicated chromosomal locations.")

            eQTLcross@model@pY <- eQTLcross@model@pY + n
            newVertexLabels <- paste0("g", seq(eQTLcross@model@pY-n+1, eQTLcross@model@pY, by=1))
            eQTLcross@model@g <- graph::addNode(newVertexLabels, eQTLcross@model@g)

            if (eQTLcross$model$p - n == 0) { ## first vertices added to the graph
              graph::nodeDataDefaults(eQTLcross@model@g, "type") <- "continuous"
              graph::nodeDataDefaults(eQTLcross@model@g, "chr") <- NA_integer_
              graph::nodeDataDefaults(eQTLcross@model@g, "loc") <- NA_real_
              graph::edgeDataDefaults(eQTLcross@model@g, "a") <- NA_real_
             }
            graph::nodeData(eQTLcross@model@g, newVertexLabels, "type") <- "continuous"
            graph::nodeData(eQTLcross@model@g, newVertexLabels, "chr") <- chr
            graph::nodeData(eQTLcross@model@g, newVertexLabels, "loc") <- location
            eQTLcross@model@vtype <- vtype <- factor(unlist(graph::nodeData(eQTLcross@model@g, nodes(eQTLcross@model@g), "type"), use.names=FALSE))
            names(eQTLcross@model@vtype) <- graph::nodes(eQTLcross@model@g)
            pY <- eQTLcross@model@pY

            eQTLcross@genes <- rbind(eQTLcross@genes, matrix(c(chr, location), nrow=n, byrow=TRUE,
                                                             dimnames=list(newVertexLabels, c("chr", "location"))))
            cVertexLabels <- graph::nodes(eQTLcross@model@g)[vtype == "continuous"]
            a <- rep(NA, pY)
            names(a) <- cVertexLabels

            eQTLcross@model@sigma <- new("dspMatrix", Dim=as.integer(c(pY, pY)),
                                         Dimnames=list(cVertexLabels, cVertexLabels),
                                         x=diag(pY)[upper.tri(diag(pY), diag=TRUE)])
            eQTLcross
          })

setMethod("addeQTL", signature(eQTLcross="eQTLcross"),
          function(eQTLcross, genes, chr=NULL, location=NULL) {
            if (eQTLcross$model$pY == 0)
              stop("Input 'eQTLcross' object has no genes to which an eQTL could be added. Use 'addGene()' first.\n")

            cVertexLabels <- graph::nodes(eQTLcross@model@g)[eQTLcross@model@vtype == "continuous"]
            mt <- match(genes, cVertexLabels)
            if (any(is.na(mt)))
              stop(sprintf("Gene(s) %s do(es) not form part of the input 'eQTLcross' object.\n", genes[is.na(mt)]))

            n <- length(genes)

            if (is.null(chr))
              chr <- sample(seq(along=qtl::chrlen(eQTLcross@map)), size=n, replace=TRUE)
            else if (length(chr) == 1)
              chr <- rep(chr, n)
            else if (length(chr) != n)
              stop("Chromosomes in 'chr' should match the number of genes.")

            if (is.null(location))
              location <- sapply(chr, function(chr) runif(1, min=0, max=qtl::chrlen(eQTLcross@map)[chr]))
            else if (length(location) != n)
              stop("Chromosomal locations in 'location' should match the number of genes.")

            ## order input eQTL to be added by genome location
            o <- order(chr, location)
            genes <- genes[o]
            chr <- chr[o]
            location <- location[o]

            ## get eQTLs that already exist
            allqtl <- alleQTL(eQTLcross)
            allqtl.loc <- paste(allqtl$chrom, allqtl$location, sep="_")
            thisqtl.loc <- paste(chr, location, sep="_")
            if (any(duplicated(paste(thisqtl.loc, genes, sep="_"))))
              stop("There are duplicated input eQTLs")
            mt.cur <- match(thisqtl.loc, allqtl.loc) ## match to current QTLs

            ## account for duplicated input QTL
            mask.newqtl <- is.na(mt.cur) & !duplicated(thisqtl.loc)

            ## create new QTLs
            n.newqtl <- sum(mask.newqtl)
            newVertexLabels <- character()
            if (n.newqtl > 0) {
              eQTLcross@model@pI <- eQTLcross@model@pI + n.newqtl
              newVertexLabels <- paste0("QTL", seq(eQTLcross@model@pI-n.newqtl+1, eQTLcross@model@pI, by=1))
              eQTLcross@model@g <- graph::addNode(newVertexLabels, eQTLcross@model@g)
              newqtlchr <- do.call("names<-", list(chr[mask.newqtl], newVertexLabels))
              newqtlloc <- do.call("names<-", list(location[mask.newqtl], newVertexLabels))
              graph::nodeData(eQTLcross@model@g, newVertexLabels, "type") <- "discrete"
              graph::nodeData(eQTLcross@model@g, newVertexLabels, "chr") <- newqtlchr[mask.newqtl]
              graph::nodeData(eQTLcross@model@g, newVertexLabels, "loc") <- newqtlloc[mask.newqtl]
              eQTLcross@model@vtype <- factor(unlist(graph::nodeData(eQTLcross@model@g, graph::nodes(eQTLcross@model@g), "type"), use.names=FALSE))
              names(eQTLcross@model@vtype) <- graph::nodes(eQTLcross@model@g)
            }

            ## find out current and new QTL from given genes
            qtlVertexLabels <- rep(NA_character_, n)
            qtlVertexLabels[!is.na(mt.cur)] <- allqtl$QTL[mt.cur[!is.na(mt.cur)]]
            qtlVertexLabels[mask.newqtl] <- newVertexLabels
            qtlVertexLabels[is.na(mt.cur) & duplicated(thisqtl.loc)] <-
              newVertexLabels[match(thisqtl.loc[is.na(mt.cur) & duplicated(thisqtl.loc)], thisqtl.loc[mask.newqtl])]

            ## add eQTLs
            eQTLcross@model@g <- graph::addEdge(from=genes, to=qtlVertexLabels, ## gi~QTLj !!
                                                graph=eQTLcross@model@g)
            graph::edgeData(eQTLcross@model@g, from=qtlVertexLabels, to=genes, "a") <- NA_real_

            eQTLcross
          })

setMethod("addGeneAssociation", signature(eQTLcross="eQTLcross"),
          function(eQTLcross, gene1, gene2) {
            if (eQTLcross$model$pY <= 1)
              stop("Input 'eQTLcross' object has less than two genes between which an association could be added. Use 'addGene()' first.\n")

            cVertexLabels <- graph::nodes(eQTLcross@model@g)[eQTLcross@model@vtype == "continuous"]
            if (is.na(match(gene1, cVertexLabels)))
              stop(sprintf("gene %s does not form part of the input 'eQTLcross' object.\n", gene1))
            if (is.na(match(gene2, cVertexLabels)))
              stop(sprintf("gene %s does not form part of the input 'eQTLcross' object.\n", gene2))
            if (gene1 == gene2)
              stop("Genes in argument 'gene1' and 'gene2' must be different.\n")

            eQTLcross@model@g <- graph::addEdge(from=gene1, to=gene2, graph=eQTLcross@model@g)

            eQTLcross
          })

## $ accessor operator
setMethod("$", signature(x="eQTLcross"),
          function(x, name) {
              switch(name,
                     type=x@type,
                     map=x@map,
                     model=x@model,
                     g=x@model$g,
                     stop("unknown eQTLcross slot or parameter. Use names() to find out which are the valid ones.")
                     )
          })

## names method
setMethod("names", signature(x="eQTLcross"),
          function(x) {
            c("g", "map", "model", "type")
          })

## alleQTL method
setMethod("alleQTL", signature(x="eQTLcross"),
          function(x) {
            eQTLedges <- graph::edges(x@model$g)[x@model$I]
            n.eqtl <- length(unlist(graph::edges(x@model$g)[x@model$I], use.names=FALSE))
            df <- as.data.frame(matrix(NA, nrow=0, ncol=5, dimnames=list(NULL, c("chrom", "location", "QTL", "gene", "a"))))

            if (n.eqtl > 0) {
              qtlvertices <- rep(names(eQTLedges), times=sapply(eQTLedges, length))
              genevertices <- unlist(eQTLedges, use.names=FALSE)
              chr <- unlist(nodeData(x@model$g, n=qtlvertices, attr="chr"), use.names=FALSE)
              loc <- unlist(nodeData(x@model$g, n=qtlvertices, attr="loc"), use.names=FALSE)
              a <- unlist(edgeData(x@model$g, from=qtlvertices, to=genevertices, attr="a"),
                          use.names=FALSE)
              df <- data.frame(chrom=chr, location=loc, QTL=qtlvertices, gene=genevertices, a=a,
                               stringsAsFactors=FALSE)
              df <- df[order(df$chrom, df$location), ]
            }

            df
          })

## ciseQTL method
setMethod("ciseQTL", signature(x="eQTLcross", cisr="missing"),
          function(x, cisr) {
            ciseQTL(x, 1)
          })

setMethod("ciseQTL", signature(x="eQTLcross", cisr="numeric"),
          function(x, cisr) {
            alleqtl <- alleQTL(x)

            if (any(is.na(match(alleqtl$gene, rownames(x@genes)))))
              stop("some eQTL genes do not match to annotated genes.")

            ciseqtl <- cbind(alleqtl, genechrom=x@genes[alleqtl$gene, 1],
                             geneloc=x@genes[alleqtl$gene, 2])
            ciseqtl <- ciseqtl[ciseqtl$chrom == ciseqtl$genechrom &
                               abs(ciseqtl$location-ciseqtl$geneloc) <= cisr, ]

            ciseqtl[, 1:5]
          })

setMethod("transeQTL", signature(x="eQTLcross", cisr="missing"),
          function(x, cisr) {
            transeQTL(x, 1)
          })

setMethod("transeQTL", signature(x="eQTLcross", cisr="numeric"),
          function(x, cisr) {
            alleqtl <- alleQTL(x)

            if (any(is.na(match(alleqtl$gene, rownames(x@genes)))))
              stop("some eQTL genes do not match to annotated genes.")

            transeqtl <- cbind(alleqtl, genechrom=x@genes[alleqtl$gene, 1], geneloc=x@genes[alleqtl$gene, 2])
            transeqtl <- transeqtl[transeqtl$chrom != transeqtl$genechrom | abs(transeqtl$location-transeqtl$geneloc) > cisr, ]
            rownames(transeqtl) <- 1:nrow(transeqtl)

            transeqtl[, 1:5]
          })

## show method
setMethod("show", signature(object="eQTLcross"),
          function(object) {
            strtype <- switch(object$type,
                              bc="backcross",
                              f2="F2",
                              NA)

            neQTLedges <- length(unlist(graph::edges(object@model$g)[object@model$I], use.names=FALSE))
            cat(sprintf("\n  eQTL %s model with %d markers, %d genes,\n  %d eQTL and %d gene-gene expression associations.\n\n",
                        strtype, qtl::totmar(object$map), object$model$pY, neQTLedges,
                        graph::numEdges(object$model$g)-neQTLedges))

            invisible(object)
          })

## plot method
setMethod("plot", signature(x="eQTLcross"),
          function(x, type=c("dot", "network"), xlab="eQTL location", ylab="Gene location", ...) {
            type <- match.arg(type)

            if (type == "network") {
              gmm <- x@model
              if (gmm$pI > 1) {
                orderedVtc <- names(sort(unlist(graph::nodeData(gmm@g, attr="loc"))[gmm@vtype == "discrete"]))
                allchr <- unlist(graph::nodeData(gmm@g, attr="chr"))[orderedVtc]
                for (chr in unique(allchr)) {
                  if (sum(allchr == chr) > 1) {
                    orderedVtcByChr <- orderedVtc[allchr == chr]
                    gmm@g <- graph::addEdge(from=orderedVtcByChr[1:(length(orderedVtc)-1)], to=orderedVtcByChr[2:length(orderedVtc)], gmm@g)
                  }
                }
              }
              plot(gmm, ...)
            } else {
              eqtl <- alleQTL(x)
              eqtlgenes <- x@genes[unique(eqtl$gene), , drop=FALSE]
              qtl <- unique(eqtl[, 1:3])
              rownames(qtl) <- qtl[, 3]
              qtl <- qtl[, 1:2]

              ## re-scale genes and marker locations in eQTL, between 0 and 1
              chrLen <- sapply(x@map, max)
              chrRelLen <- chrLen / sum(chrLen)
              chrRelCumLen <- c(0, cumsum(chrRelLen)[-length(chrRelLen)])
              geneRelPos <- chrRelCumLen[eqtlgenes[, "chr"]] +
                            chrRelLen[eqtlgenes[, "chr"]]*(eqtlgenes[, "location"] / chrLen[eqtlgenes[, "chr"]])
              names(geneRelPos) <- rownames(eqtlgenes)

              qtlRelPos <- chrRelCumLen[qtl[, "chrom"]] +
                           chrRelLen[qtl[, "chrom"]]*(qtl[, "location"]/chrLen[qtl[, "chrom"]])
              names(qtlRelPos) <- rownames(qtl)

              plot(qtlRelPos[eqtl[, "QTL"]], geneRelPos[eqtl[, "gene"]], pch=".",
                   xlab=xlab, ylab=ylab, cex=4, axes=FALSE, ...)
              segments(c(chrRelCumLen, 1), 0, c(chrRelCumLen, 1), 1, col=gray(0.75), lty="dotted", lwd=2)
              segments(0, c(chrRelCumLen, 1), 1, c(chrRelCumLen, 1), col=gray(0.75), lty="dotted", lwd=2)
              axis(1, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
                   labels=names(chrLen), tick=FALSE, cex.axis=0.9)
              axis(2, at=(chrRelCumLen + c(chrRelCumLen[-1], 1))/2,
                   labels=names(chrLen), tick=FALSE, cex.axis=0.9, las=1)
            }
          }) 

## eQTLcrossParam constructor
eQTLcrossParam <- function(map=do.call("class<-", list(list("1"=do.call("class<-", list(do.call("names<-", list(seq(0, 100, length.out=20), paste0("m", 1:20))), "A"))), "map")),
                           type="bc", genes=20, cis=1, trans=as.integer(NULL), cisr=1, d2m=0,
                           networkParam=dRegularGraphParam()) {
  type <- match.arg(type)

  if (class(map) == "list")
    class(map) <- "map"

  if (!is(networkParam, "graphParam"))
    stop("argument 'networkParam' should be a 'graphParam' object defining the type of gene network to simulate.")

  ## if the number of continuous variables in networkParam does not match the number
  ## of genes then the latter overrides the number in networkParam
  nm <- qtl::totmar(map)
  if (networkParam@p != genes)
    networkParam@p <- as.integer(genes)
  networkParam@labels <- paste0("g", 1:genes)

  if (cis < 0)
    stop("argument 'cis' should be either a real number between 0 and 1 representing the wanted fraction of genes with cis-eQTL associations or a positive integer number specifying the number of desired cis-eQTL associations.")

  n.cisQTL <- ifelse(cis > 1, cis, floor(networkParam@p * cis)) ## number of cisQTL
  n.transQTL <- length(trans)                                   ## number of transQTL

  if (n.cisQTL+n.transQTL > nm && d2m == 0)
    stop(sprintf("more eQTL (%d) than markers (%d). Either, increase marker density in the genetic map, increase d2m or decrease the number of eQTL.", n.cisQTL+n.transQTL, nm))

  if (n.cisQTL > nm && d2m == 0)
    stop("more cis-QTL than markers. Either, increase marker density in the genetic map, increase d2m or decrease the number of cis-QTL.")

  new("eQTLcrossParam", map=map, type=type, cis=cis, trans=as.integer(trans), cisr=cisr, d2m=d2m,
      networkParam=networkParam)
}

## eQTLcrossParam show method
setMethod("show", signature(object="eQTLcrossParam"),
          function(object) {
            strtype <- switch(object@type,
                              bc="backcross",
                              f2="F2",
                              NA)

            n.cisQTL <- ifelse(object@cis > 1, object@cis, floor(object@cis*object@networkParam@p))
            cat(sprintf("\n  eQTL %s parameter object defining %d markers,\n",
                        strtype, qtl::totmar(object@map)))
            cat(sprintf("  %d genes, %d cis-eQTL and %d trans-eQTL.\n\n", object@networkParam@p,
                        n.cisQTL, sum(object@trans)))
            cat(sprintf("  cis-eQTL associations occur within %.1f cM of a gene\n", object@cisr))
            cat(sprintf("  and all eQTL are located at %.1f cM from a marker.\n\n", object@d2m))
            cat("  Gene network parameters are defined by a\n")
            print(object@networkParam)

            invisible(object)
          })

## constructor simulation methods
## genes: # of genes or vector of cM positions, each element corresponding to distinct gene
## cis: fraction of genes with eQTL in cis
## trans: vector of integer numbers, each element corresponding to a trans-eQTL and its value
##        is the number of genes linked to this eQTL
## cisr: cis radius, maximum distance in cM that defines a QTL in cis wrt to a gene
## d2m: distance of every gene or eQTL to a marker, default 0 implies eQTLs and genes are located at markers,
##      and therefore, the maximum number of eQTL is bounded by the number of markers
setMethod("reQTLcross", signature(n="eQTLcrossParam", network="missing"),
          function(n, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=1L, network=n, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="missing", network="eQTLcrossParam"),
          function(n, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=1L, network, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="numeric", network="eQTLcrossParam"),
          function(n=1, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=as.integer(n), network, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="integer", network="eQTLcrossParam"),
          function(n=1L, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            map <- network@map
            nm <- qtl::totmar(map)                 ## total number of markers
            nmbychr <- qtl::nmar(map)              ## number of markers per chromosome
            csnmbychr <- cumsum(nmbychr)           ## cumulative sum of the number of markers per chromosome
            mcmloc <- unlist(map, use.names=FALSE) ## location of markers in cM
            cmlenbychr <- sapply(map, max)         ## chromosome length in cM
            cscmlenbychr <- cumsum(cmlenbychr)     ## cumulative sum of chromosome length in cM

            type <- network@type
            pY <- network@networkParam@p
            Y <- network@networkParam@labels

            cis <- network@cis
            trans <- network@trans
            d2m <- network@d2m
            cisr <- network@cisr

            idx.m.cisQTL <- chr.genes.cisQTL <- loc.genes.cisQTL <- c()
            ## simulate gene locations in cM for genes with cis-QTL

            n.genes <- pY
            n.cisQTL <- ifelse(cis > 1, cis, floor(n.genes * cis)) ## number of cisQTL
            n.transQTL <- length(trans)                            ## number of transQTL

            if ((class(a) == "numeric" || class(a) == "integer") && length(a) > 1 && length(a) != n.cisQTL + n.transQTL)
              stop(sprintf("argument 'a' contains %d values of eQTL additive effects while arguments 'genes', 'cis' and 'trans' determine a total number of %d eQTL.", length(a), n.cisQTL+n.transQTL))

            if (class(a) == "function" && length(formals(a)) != 1)
              stop("when argument 'a' is a function it should contain one argument taking the number of eQTL.")

            sim <- list()
            for (i in 1:n) {

              chr.genes.cisQTL <- loc.genes.cisQTL <- numeric(0)
              chr.genes.nocisQTL <- loc.genes.nocisQTL <- numeric(0)

              if (n.cisQTL > 0) { ## place genes with cis-eQTL associations
                nocis <- FALSE
                j <- 0
                ## enforce genes with cis-QTL being located at least 1 x cisr cM apart
                while (!nocis && j < 10) {
                  idx.m.cisQTL <- sample(1:nm, size=n.cisQTL, replace=FALSE)
                  chr.genes.cisQTL <- sapply(idx.m.cisQTL, function(j, cs) sum(cs < j)+1, csnmbychr)
                  loc.genes.cisQTL <- mcmloc[idx.m.cisQTL] + d2m
                  nocis <- sapply(split(loc.genes.cisQTL, chr.genes.cisQTL),
                                  function(gxc, cisr) {
                                    nocis <- TRUE
                                    if (length(gxc) > 1)
                                      nocis <- all(combn(gxc, 2, function(x) abs(x[1]-x[2])) > 1*cisr)
                                    nocis
                                    }, cisr)
                  nocis <- all(nocis)
                  j <- j + 1
                }

                if (!nocis)
                  stop("impossible to simulate genes with cis-eQTL. Either decrease cisr, decrease the number of genes, or increase marker density in the genetic map.")
              }

              if (n.genes - n.cisQTL > 0) { ## place genes with trans-eQTL associations
                n.genes.left <- n.genes - n.cisQTL

                genes.cisQTL <- do.call("names<-", list(vector(mode="list", length=length(nmbychr)),
                                                        names(nmbychr))) ## just in case n.cisQTL == 0
                if (n.cisQTL > 0)
                  genes.cisQTL <- split(loc.genes.cisQTL, chr.genes.cisQTL)
                chr.genes.nocisQTL <- sample(1:length(nmbychr), size=n.genes.left, replace=TRUE)
                loc.genes.nocisQTL <- sapply(1:n.genes.left,
                                             function(j, chr, chrlen, genes.cisQTL, cisr) {
                                               k <- 0
                                               nocis <- FALSE
                                               pos <- NA
                                               while (!nocis && k < 10) {
                                                 pos <- runif(1, min=0, max=chrlen[chr[j]])
                                                 if (!is.na(match(as.character(chr[j]), names(genes.cisQTL))))
                                                   nocis <- all(abs(pos-genes.cisQTL[[as.character(chr[j])]]) > 1*cisr)
                                                 else
                                                   nocis <- TRUE
                                                 k <- k + 1
                                               }
                                               if (!nocis)
                                                 stop(sprintf("impossible to simulate %d genes without cis-eQTL and located %.1fcM away from the %d cis-eQTL genes. Either decrease cisr, decrease the total number of genes or those with cis-eQTL.", n.genes.left, cisr, n.cisQTL))

                                               pos
                                             }, chr.genes.nocisQTL, cmlenbychr, genes.cisQTL, cisr)
              }

              ## build gene annotation matrix
              genes <- rbind(cbind(chr.genes.cisQTL, loc.genes.cisQTL),
                             cbind(chr.genes.nocisQTL, loc.genes.nocisQTL))
              o <- order(genes[, 1], genes[, 2])
              colnames(genes) <- c("chr", "location")
              rownames(genes)[o] <- Y ## name genes by their location along the genome
              cisQTLgenes <- character(0)
              if (n.cisQTL > 0)
                cisQTLgenes <- rownames(genes)[1:n.cisQTL] ## take first their gene ID
              genes <- genes[o, , drop=FALSE] ## reorder genes by their location along the genome
              nocisQTLgenes <- setdiff(Y, cisQTLgenes)
              cisQTLgenes <- match(cisQTLgenes, rownames(genes))
              nocisQTLgenes <- match(nocisQTLgenes, rownames(genes))
              cisQTL <- matrix(NA_real_, nrow=0, ncol=3)
              if (n.cisQTL > 0)
                cisQTL <- cbind(genes[cisQTLgenes, "chr"], genes[cisQTLgenes, "location"], cisQTLgenes)

              ## simulate trans-QTL associations

              ## function to search markers (m) in cis to a gene (g) within a radius (r)
              cism <- function(markers, gene, radius) which(markers >= gene-radius & markers <= gene+radius)

              transQTL <- matrix(NA, nrow=0, ncol=3)
              if (sum(trans) > length(nocisQTLgenes))
                stop(sprintf("not enough genes (%d) to simulate %d trans-QTL associations. Either decrease the number of trans-QTL associations or decrease 'cis'.", length(nocisQTLgenes), sum(trans)))

              for (ng in trans) {
                tmpmap <- map
                transQTLgenes <- sample(nocisQTLgenes, size=ng, replace=FALSE)
                nocisQTLgenes <- setdiff(nocisQTLgenes, transQTLgenes) ## removed sampled trans-genes
                transQTLgenes <- split(transQTLgenes, names(tmpmap)[genes[transQTLgenes, "chr"]])
                for (chr in names(transQTLgenes)) {
                  loc.genes <- genes[transQTLgenes[[chr]], "location"]
                  allcm <- unlist(sapply(loc.genes,
                                         function(g, m, cisr) cism(m, g, cisr),
                                           tmpmap[[chr]]+d2m, cisr)) ## all cis-markers
                  if (length(allcm) > 0)
                    tmpmap[[chr]] <- tmpmap[[chr]][-allcm] ## remove all cis-markers, if found
                }
                tmpmap2 <- tmpmap
                for (chr in 1:length(tmpmap)) {
                     tmpmap2[[chr]] <- tmpmap2[[chr]][is.na(match(tmpmap[[chr]], transQTL[transQTL[, 1] == chr, 2]))]
                     class(tmpmap2[[chr]]) <- class(tmpmap[[chr]])
                }
                tmpmap <- tmpmap2

                tmpmcmloc <- unlist(tmpmap, use.names=FALSE) + d2m
                tmpnmbychr <- qtl::nmar(tmpmap)
                tmpcsnmbychr <- cumsum(tmpnmbychr)
                if (sum(tmpnmbychr) < 1)
                  stop("no markers left to sample a marker hotspot. Either decrease the number of trans-QTL associations or decrease 'cisr'.")
                tm <- sample(sum(tmpnmbychr), size=1, replace=FALSE) ## sample the hotspot marker
                chr.tm <- sum(tmpcsnmbychr < tm) + 1
                transQTL <- rbind(transQTL, cbind(rep(chr.tm, times=ng), rep(tmpmcmloc[tm], times=ng),
                                                  unlist(transQTLgenes, use.names=FALSE)))
              }
    
              ## simulate gene network
              sim.g <- rgraphBAM(n=1L, network@networkParam, verbose=verbose)
              edg <- graph::edges(sim.g)
              edg <- cbind(rep(names(edg), times=sapply(edg, length)),
                          unlist(edg, use.names=FALSE))
              if (nrow(edg) > 0)
                edg <- unique(t(apply(edg, 1, sort)))
              edg <- cbind(match(edg[, 1], Y), match(edg[, 2], Y))

              ## simulate conditional covariance matrix
              sim.sigma <- qpG2Sigma(g=sim.g, rho=rho, tol=tol, verbose=verbose)

              ## build matrix of eQTL associations, genes are indexed wrt to the rownames of the 'genes' matrix
              qtl <- rbind(cisQTL, transQTL)

              ## simulate additive effect in QTL
              if ((class(a) == "numeric" || class(a) == "integer") && length(a) == 1)
                a <- rep(a, times=nrow(qtl))
              else if (class(a) == "function")
                a <- a(nrow(qtl))

              qtl <- cbind(qtl, a)

              sim[[i]] <- eQTLcross(map, genes=genes, model=qtl, type=type,
                          geneNetwork=edg, rho=rho, sigma=sim.sigma)
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

setMethod("reQTLcross", signature(n="eQTLcross", network="missing"),
          function(n, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=1L, network=n, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="missing", network="eQTLcross"),
          function(n, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=1L, network, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="numeric", network="eQTLcross"),
          function(n=1, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {
            reQTLcross(n=as.integer(n), network, rho, a, tol, verbose)
          })

setMethod("reQTLcross", signature(n="integer", network="eQTLcross"),
          function(n=1L, network, rho=0.5, a=1, tol=0.001, verbose=FALSE) {

            ## simulated gene network is fixed by the input argument 'network'
            genes <- network$model$Y
            eQTLs <- graph::edges(network$g)[network$model$I]
            n.eQTLs <- sum(sapply(eQTLs, length))
            sim.g <- graph::subGraph(genes, network$g)
            Ylabels <- network$model$Y ## to be used later to re-order sigma rows and columns

            if ((class(a) == "numeric" || class(a) == "integer") && length(a) > 1 && length(a) != n.eQTL)
              stop(sprintf("argument 'a' contains %d values of eQTL additive effects while the total number of eQTL is %d.", length(a), n.eQTL))

            if (class(a) == "function" && length(formals(a)) != 1)
              stop("when argument 'a' is a function it should contain one argument taking the number of eQTL.")

            if ((class(a) == "numeric" || class(a) == "integer") && length(a) == 1)
              a <- rep(a, times=n.eQTLs)
            else if (class(a) == "function")
              a <- a(n.eQTLs)

            ## additive effects per gene are the sum of additive effects from each eQTL (not yet useful!)
            ## while additive effects are stored by eQTL it is still not possible to specify different
            ## effects to different eQTLs
            network@model@a <- do.call("names<-", list(rep(0, length(genes)), genes))
            k <- 1
            for (i in seq(along=eQTLs)) {
              graph::edgeData(network@model@g, from=eQTLs[[i]],
                              to=rep(names(eQTLs)[i], length(eQTLs[[i]])), "a") <- a[k]
              network@model@a[eQTLs[[i]]] <- network@model@a[eQTLs[[i]]] + a[k]
              k <- k + 1
            }

            network@model@rho <- rho

            sim <- list()
            for (i in 1:n) {
    
              ## simulate covariance matrix to be interpreted conditionally later
              sim.sigma <- qpG2Sigma(g=sim.g, rho=rho, tol=tol, verbose=verbose)
              rownames(sim.sigma) <- colnames(sim.sigma) <- nodes(sim.g)
              sim.sigma <- sim.sigma[Ylabels, Ylabels, drop=FALSE] ## put back rows and columns into the original variable order
                                                                   ## since 'subgraph()' re-orders nodes alphabetically (sigh!)

              network@model@sigma <- sim.sigma

              sim[[i]] <- network
            }

            if (n == 1)
              sim <- sim[[1]]

            sim
          })

## overload sim.cross() from the qtl package to enable simulating eQTL data from
## an experimental cross

sim.cross <- function(map, model, ...) UseMethod("sim.cross", model)
sim.cross.default <- function(map, model, ...) qtl::sim.cross(map, model, ...)
sim.cross.matrix <- function(map, model, ...) qtl::sim.cross(map, model, ...)
setMethod("sim.cross", c(map="map", model="matrix"), sim.cross.matrix)

sim.cross.eQTLcross <- function(map, model, n.ind=100, ...) {

            eQTLs <- alleQTL(model)[, c("chrom", "location")]
            crossModel <- unique(eQTLs)

            if (is.null(model@model@a) ||
                (all(model@model$sigma[upper.tri(model@model$sigma)] == 0) && graph::numEdges(model$g) - nrow(eQTLs) > 0))
              stop("Parameters for the input 'eQTLcross' model have not been simulated. Please simulated them running 'reQTLcross()' with this 'eQTLcross' model as input argument.")

            crossModel <- cbind(crossModel, dummya=rep(1, nrow(crossModel)))
            cross <- qtl::sim.cross(map=map, model=crossModel, type=model$type, n.ind=n.ind, ...)

            ## we use the same data matrix X to store the mean values employed
            ## during the simulation process.
            I <- model@model$I
            Y <- model@model$Y
            stopifnot(identical(I, colnames(cross$qtlgeno)))
            cross$pheno <- qpgraph:::calculateCondMean(model@model, cross$qtlgeno) 
            rownames(cross$pheno) <- 1:n.ind
            colnames(cross$pheno) <- Y
            cross$pheno <- as.data.frame(cross$pheno)

            YxI <- Y[which(sapply(graph::edges(model@model$g)[Y], function(xYi, vt) sum(vt[xYi] == "discrete"), model@model@vtype) > 0)]
            xtab <- tapply(1:n.ind, apply(cross$pheno[, YxI, drop=FALSE], 1, function(i) paste(i, collapse="")))
            xtab <- split(as.data.frame(cross$qtlgeno[, I]), xtab)
            for (i in 1:length(xtab)) {
              li <- xtab[[i]]
              which_n <- as.numeric(rownames(li))
              cross$pheno[which_n, Y] <- mvtnorm::rmvnorm(length(which_n), mean=as.numeric(cross$pheno[which_n[1], Y]),
                                                          sigma=as.matrix(model@model$sigma))
            }

            cross
          }
setMethod("sim.cross", signature(map="map", model="eQTLcross"), sim.cross.eQTLcross)
