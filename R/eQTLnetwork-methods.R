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
            cat(sprintf("  Input size: %d markers %d genes\n", length(mNames), length(gNames)))
            cat("  Model formula: ", deparse(object@modelFormula), "\n")
            g <- object@qpg@g
            if (numNodes(g) > 0) {
              edg <- extractFromTo(g)
              qstr <- NULL
              if (!is.na(object@p.value))
                qstr <- "0"
              if (!is.na(object@epsilon))
                qstr <- ifelse(qstr == "0",
                               sprintf("0,%s", paste(object@qOrders, collapsed=",")),
                               paste(object@qOrders, collapsed=","))
              if (!is.na(object@alpha))
                qstr <- paste(qstr, "*", sep=",")
              cat(sprintf("  G^(%s): %d vertices and %d edges corresponding to\n         %d eQTLs and %d gene-gene associations\n         meeting a %s-adjusted p-value < %.2f\n         and involving %d genes and %d loci\n",
                          qstr, numNodes(g), numEdges(g),
                          sum(edg$from %in% mNames | edg$to %in% mNames), sum(edg$from %in% gNames & edg$to %in% gNames),
                          object@adjustMethod, object@p.value,
                          length(intersect(nodes(g), gNames)), length(intersect(nodes(g), mNames))))
            }
          })
