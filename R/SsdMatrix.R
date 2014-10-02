setValidity("SsdMatrix",
            function(object) {
              x <- prod(dim(object@ssd))
              if (!is.na(object@n) && x > 0 && object@n < 1)
                "'n' should equal the number of observations of the data employed to estimate the ssd"
              else if (!is.na(object@n) && x == 0 && object@n > 0)
                "'n' should equal 0 when ssd is empty"
              else
                TRUE
            })

## adapted from Matrix:::prMatrix()
setMethod("show", signature(object="SsdMatrix"),
          function(object) {
            digits <- getOption("digits")
            maxp <-getOption("max.print")
            d <- dim(object)
            cl <- class(object)
            tri <- extends(cl, "triangularMatrix")
            xtra <- if(tri && object@diag == "U") " (unitriangular)" else ""
            cat(sprintf('%d x %d SsdMatrix extending class "%s"%s\n',
                        d[1], d[2], cl, xtra))
            if (prod(d) <= maxp) {
              if (tri) {
                upper <- object@uplo == "U"
                m <- as(object, "matrix")
                cf <- format(m, digits=digits, justify="none")
                cf[if (upper) row(cf) > col(df) else row(cf) < col(cf)] <- "."
                print(cf, quote=FALSE, right=TRUE, max=maxp)
              } else
                print(as(object, "matrix"), digits=digits, max=maxp)
            } else { ## d[1] > maxp / d[2] >= nr
              m <- as(object, "matrix")
              nr <- maxp %/% d[2]
              n2 <- ceiling(nr / 2)
              print(head(m, max(1, n2)))
              cat("\n...........\n\n")
              print(tail(m, max(1, nr - n2)))
              cat("\n...........\n\n")
            }
            if (!is.na(object@n))
              cat(sprintf("n = %d\n", object@n))
            else
              cat("n = NA\n")

            invisible(object)
          })

setMethod("determinant", signature(x = "SsdMatrix", logarithm = "missing"),
          function(x, logarithm, ...)
            Matrix::determinant(x@ssd, logarithm=TRUE, ...))

setMethod("dim", signature(x = "SsdMatrix"),
          function(x)
            x@ssd@Dim, valueClass = "integer")

setMethod("dimnames", signature(x = "SsdMatrix"),
          function(x)
            x@ssd@Dimnames)

setAs("SsdMatrix", "matrix",
      function(from)
        as(from@ssd, "matrix"))

setAs("SsdMatrix", "dspMatrix",
      function(from)
        from@ssd)
