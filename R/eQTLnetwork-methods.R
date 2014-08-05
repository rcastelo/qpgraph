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
            cat(sprintf("  %d markers %d genes\n",
                max(c(length(unlist(object@geneticMap, use.names=FALSE)),
                      length(unlist(object@physicalMap, use.names=FALSE)))),
                length(object@geneAnnotation)))
            cat("  Model formula: ", deparse(object@modelFormula), "\n")
          })
