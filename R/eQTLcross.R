setValidity("eQTLcross",
            function(object) {
              valid <- TRUE

              if (object@type == "bc" && object@model@dLevels != 2)
                valid <- "when 'type'=\"bc\" the 'HMgmm' object should have 2 levels for the discrete variables."

              valid
            })

## constructor methods
setMethod("eQTLcross", signature(map="map", genes="missing", model="missing"),
          function(map, genes, model, type="bc",
                   nGenes=100, geneNetwork=matrix(NA, nrow=0, ncol=2),
                   rho=0.5, sigma=diag(nGenes)) {
            genes <- matrix(NA, nrow=0, ncol=2)
            model <- matrix(NA, nrow=0, ncol=4)
            eQTLcross(map, genes, model, type, nGenes, geneNetwork, rho, sigma)
          })

## model is a matrix where each row corresponds to a different eQTL, and
## gives the chromosome number, cM position, gene and additive effect of the eQTL
setMethod("eQTLcross", signature(map="map", genes="missing", model="matrix"),
          function(map, genes, model, type="bc",
                   nGenes=100, geneNetwork=matrix(NA, nrow=0, ncol=2),
                   rho=0.5, sigma=diag(nGenes)) {
            genes <- matrix(NA, nrow=0, ncol=2)
            eQTLcross(map, genes, model, type, nGenes, geneNetwork, rho, sigma)
          })

setMethod("eQTLcross", signature(map="map", genes="matrix", model="matrix"),
          function(map, genes, model, type="bc",
                   nGenes=100, geneNetwork=matrix(NA, nrow=0, ncol=2),
                   rho=0.5, sigma=diag(nGenes)) {

            if (nGenes < 1)
              stop("nGenes < 1.")

            dVertexLabels <- c()
            cVertexLabels <- paste0("g", 1:nGenes)
            edges <- as.data.frame(matrix(NA, nrow=0, ncol=3,
                                          dimnames=list(NULL, c("from", "to", "weight"))),
                                   stringsAsFactors=FALSE)
            if (nrow(model) > 0) {
              model <- model[order(model[, 1], model[, 2]), ]             ## sort eQTL by chromosome and location
              rownames(model) <- NULL                                     ## remove rownames, if any
              eQTLtag <- rle(apply(model[, 1:2], 1, paste, collapse="_")) ## tag every different eQTL location
              nQTL <- length(eQTLtag$values)                              ## count the number of different QTL
              model <- cbind(model, rep(1:nQTL, times=eQTLtag$lengths))   ## add QTL number to each different QTL
              dVertexLabels <- paste0("QTL", 1:nQTL)                      ## assign QTL names to each different QTL
              edges <- data.frame(from=paste0("QTL", model[, 5]),
                                  to=paste0("g", model[, 3]),
                                  weight=1, stringsAsFactors=FALSE)
            }
            
            if (nrow(geneNetwork) > 0)
              edges <- rbind(edges,
                             data.frame(from=paste0("g", geneNetwork[, 1]),
                                        to=paste0("g", geneNetwork[, 2]),
                                        weight=rep(1, nrow(geneNetwork)),
                                        stringsAsFactors=FALSE))

            g <- graph::graphBAM(edges, nodes=c(dVertexLabels, cVertexLabels))
            graph::nodeDataDefaults(g, "type") <- "continuous"
            graph::nodeDataDefaults(g, "chr") <- NA_integer_
            graph::nodeDataDefaults(g, "loc") <- NA_real_
            graph::edgeDataDefaults(g, "a") <- NA_real_

            if (nrow(model) > 0) {
              qtlchr <- unique(model[, c(1, 5)])
              qtlloc <- unique(model[, c(2, 5)])
              qtlchr <- do.call("names<-", list(qtlchr[, 1], paste0("QTL", qtlchr[, 2])))
              qtlloc <- do.call("names<-", list(qtlloc[, 1], paste0("QTL", qtlloc[, 2])))
              graph::nodeData(g, dVertexLabels, "type") <- "discrete"
              graph::nodeData(g, dVertexLabels, "chr") <- qtlchr[dVertexLabels]
              graph::nodeData(g, dVertexLabels, "loc") <- qtlloc[dVertexLabels]
              if (nrow(model) > 0) {
                mixedEdgeMask <- !is.na(match(edges$from, dVertexLabels))
                graph::edgeData(g, from=edges[mixedEdgeMask, ]$from, to=edges[mixedEdgeMask, ]$to, "a") <- model[, 4]
                ## 21/2/2013 because of the attribute bug in graphBAM by now we store additive effects in the weight attribute
                graph::edgeData(g, from=edges[mixedEdgeMask, ]$from, to=edges[mixedEdgeMask, ]$to, "weight") <- model[, 4]
              }
            }

            eQTLcross(map=map, genes=genes, model=g, type=type, rho=rho, sigma=sigma)
          })

setMethod("eQTLcross", signature(map="map", genes="matrix", model="graphBAM"),
          function(map, genes, model,
                   type="bc",
                   rho=0.5,
                   sigma=diag(sum(unlist(graph::nodeData(model, graph::nodes(model), attr="type"), use.names=FALSE) == "continuous"))) {
            type <- match.arg(type)

            dLevels <- switch(type,
                              bc=2L)

            maskcvertex <- unlist(graph::nodeData(model, graph::nodes(model), attr="type"), use.names=FALSE) == "continuous"
            cVertexLabels <- graph::nodes(model)[maskcvertex]

            ## additive effects per gene are the sum of additive effects from each eQTL
            a <- sapply(cVertexLabels, function(v) sum(unlist(edgeData(model, to=v, attr="a"), use.names=FALSE), na.rm=TRUE))

            if (is.null(rownames(sigma)) || is.null(colnames(sigma)))
              rownames(sigma) <- colnames(sigma) <- cVertexLabels

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
            qtlvertices <- rep(names(eQTLedges), times=sapply(eQTLedges, length))
            genevertices <- unlist(eQTLedges, use.names=FALSE)
            chr <- unlist(nodeData(x@model$g, n=qtlvertices, attr="chr"), use.names=FALSE)
            loc <- unlist(nodeData(x@model$g, n=qtlvertices, attr="loc"), use.names=FALSE)
            a <- unlist(edgeData(x@model$g, from=qtlvertices, to=genevertices, attr="a"), use.names=FALSE)
            ## 21/2/2013 because of the attribute bug in graphBAM by now we store additive effects in the weight attribute
            a <- unlist(edgeData(x@model$g, from=qtlvertices, to=genevertices, attr="weight"), use.names=FALSE)
            df <- data.frame(chrom=chr, location=loc, QTL=qtlvertices, gene=genevertices, a=a,
                             stringsAsFactors=FALSE)
            df <- df[order(df$chrom, df$location), ]
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

            ciseqtl <- cbind(alleqtl, genechrom=x@genes[alleqtl$gene, 1], geneloc=x@genes[alleqtl$gene, 2])
            ciseqtl <- ciseqtl[ciseqtl$chrom == ciseqtl$genechrom & abs(ciseqtl$location-ciseqtl$geneloc) <= cisr, ]

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
            transeqtl <- transeqtl[transeqtl$chrom != transeqtl$genechrom || abs(transeqtl$location-transeqtl$geneloc) > cisr, ]

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

            if (type == "network")
              plot(x@model, ...)
            else {
              eqtl <- alleQTL(x)
              eqtlgenes <- x@genes[unique(eqtl$gene), ]
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
