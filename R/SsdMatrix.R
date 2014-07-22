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

setMethod("show", signature(object = "SsdMatrix"),
          function(object) {
            Matrix:::prMatrix(object@ssd)
            if (!is.na(object@n))
              cat(sprintf("n = %d\n", object@n))
            else
              cat("n = NA\n")
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
