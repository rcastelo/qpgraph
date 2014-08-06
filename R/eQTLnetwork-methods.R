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
              mNames <- as.vector(sapply(object@geneticMap, names))
            else if (length(object@physicalMap) > 0)
              mNames <- as.vector(sapply(object@physicalMap, names))
            else
              warning("Both the genetic and the physical map are empty.")

            gNames <- names(object@geneAnnotation)
            cat(sprintf("  Total #: %d markers %d genes\n", length(mNames), length(gNames)))
            cat("  Model formula: ", deparse(object@modelFormula), "\n")
            g.0 <- object@g.0
            if (numNodes(g.0) > 0) {
              edg <- extractFromTo(g.0)
              cat(sprintf("  G^(0): %d vertices and %d edges corresponding to\n         %d eQTLs and %d gene-gene associations\n         involving %d genes and %d loci\n",
                          numNodes(g.0), numEdges(g.0),
                          sum(edg$from %in% mNames | edg$to %in% mNames), sum(edg$from %in% gNames & edg$to %in% gNames),
                          length(intersect(nodes(g.0), gNames)), length(intersect(nodes(g.0), mNames))))
            }
            g.q <- object@g.q@g
            if (numNodes(g.q) > 0) {
              edg <- extractFromTo(g.q)
              cat(sprintf("  G^(q): %d vertices and %d edges corresponding to\n         %d eQTLs and %d gene-gene associations\n         involving %d genes and %d loci\n",
                          numNodes(g.q), numEdges(g.q),
                          sum(edg$from %in% mNames | edg$to %in% mNames), sum(edg$from %in% gNames & edg$to %in% gNames),
                          length(intersect(nodes(g.q), gNames)), length(intersect(nodes(g.q), mNames))))
            }
            if (numNodes(g.0) > 0 && numNodes(g.q) > 0) {
              g <- graphIntersect(g.0, g.q)
              edg <- extractFromTo(g)
              cat(sprintf("  G^(0,q): %d vertices and %d edges corresponding to\n           %d eQTLs and %d gene-gene associations\n           involving %d genes and %d loci\n",
                          numNodes(g), numEdges(g),
                          sum(edg$from %in% mNames | edg$to %in% mNames), sum(edg$from %in% gNames & edg$to %in% gNames),
                          length(intersect(nodes(g), gNames)), length(intersect(nodes(g), mNames))))
            }
          })
