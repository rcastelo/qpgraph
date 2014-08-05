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
            mNames <- unique(c(names(unlist(object@geneticMap, use.names=FALSE)),
                               names(unlist(object@physicalMap, use.names=FALSE))))
            gNames <- names(object@geneAnnotation)
            cat(sprintf("  %d markers %d genes\n", length(mNames), length(gNames)))
            cat("  Model formula: ", deparse(object@modelFormula), "\n")
            g.0 <- object@g.0
            if (numNodes(g.0) > 0)
              cat(sprintf("  G^(0): %d eQTLs and %d gene-gene associations\n",
                          length(unlist(edges(g.0)[mNames], use.names=FALSE)),
                          length(unlist(lapply(edges(g.0)[gNames], function(x, m) intersect(x, m), gNames), use.names=FALSE))))
            g.q <- object@g.q
            if (numNodes(g.q) > 0)
              cat(sprintf("  G^(q): %d eQTLs and %d gene-gene associations\n",
                          length(unlist(edges(g.q)[mNames], use.names=FALSE)),
                          length(unlist(lapply(edges(g.q)[gNames], function(x, m) intersect(x, m), gNames), use.names=FALSE))))
            if (numNodes(g.0) > 0 && numNodes(g.q) > 0) {
              g <- graphIntersect(g.0, g.q)
              cat(sprintf("  G^(0,q): %d eQTLs and %d gene-gene associations\n",
                          length(unlist(edges(g)[mNames], use.names=FALSE)),
                          length(unlist(lapply(edges(g)[gNames], function(x, m) intersect(x, m), gNames), use.names=FALSE))))
            }
          })
