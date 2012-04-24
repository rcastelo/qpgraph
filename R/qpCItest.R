## function: qpCItest
## purpose: perform a conditional independence test between two variables given
##          a conditioning set
## parameters: X - data where to perform the test
##             i - index or name of one of the two variables in X to test
##             j - index or name of the other variable in X to test
##             Q - indexes or names of the variables in X forming the conditioning set
##             I - indexes or names of the variables in X that are discrete
##             n - sample size (when data is directly the sample covariance matrix)
##             long.dim.are.variables - if TRUE it assumes that when the data is
##                                      a data frame or a matrix, the longer
##                                      dimension is the one defining the random
##                                      variables, if FALSE then random variables
##                                      are assumed to be at the columns
##             exact.test - employ an exact test when working with HMGMs
##             R.code.only - flag set to FALSE when using the C implementation
## return: a list with two members, the t-statistic value and the p-value
##         on rejecting the null hypothesis of independence

setGeneric("qpCItest", function(X, ...) standardGeneric("qpCItest"))

## ADD gLevels to either directly denote the number of genotype levels or trigger their
## automatic calculation

## X comes as an smlSet object
setMethod("qpCItest", signature(X="smlSet"),
          function(X, i=1, j=2, Q=c(), exact.test=TRUE, R.code.only=FALSE) {
            p <- as.integer(nrow(X))
            h <- as.integer(ncol(Biobase::pData(X)))
            Xsml <- GGBase::smList(X)
            sByChr <- sapply(Xsml, ncol)
            gLevels <- sum(unique(as.vector(as(Xsml[[1]][, 1:min(sByChr[1], 1000)], "matrix"))) > 0)
            cumsum_sByChr <- c(0, cumsum(sByChr))
            s <- sum(sByChr)
            n <- as.integer(ncol(X))
            fNames <- Biobase::featureNames(X)
            pNames <- colnames(Biobase::pData(X))
            sNamesByChr <- lapply(Xsml, colnames)
            sNames <- unlist(sNamesByChr, use.names=FALSE)
            Xsub <- matrix(0, nrow=n, ncol=length(c(i, j, Q)))
            colnames(Xsub) <- 1:length(c(i, j, Q))
            nLevels <- rep(NA, length(c(i, j, Q)))
            missingMask <- rep(FALSE, length(c(i, j, Q)))
            x <- Y <- I <- c()

            nam_i <- i
            if (is.character(i)) {
              smt <- match(i, sNames)
              if (is.na(match(i, fNames)) && is.na(match(i, pNames)) && is.na(smt))
                stop(sprintf("i=%s does not form part of the variable names of the data\n", i))
              if (!is.na(match(i, fNames))) ## then 'i' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[i, ]
              else if (!is.na(match(i, pNames))) ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i]
              else { ## then 'i' refers to a SNP genotype (discrete)
                x <- as(Xsml[[sum(cumsum_sByChr < smt)]][, i], "numeric")+1
                if (any(x > 3, na.rm=TRUE))
                  warning(sprintf("i=%s has uncertain genotype calls which are treated here as missing", nam_i))
                x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
                nLevels[1] <- gLevels
              }
            } else {
              if (i <= p) ## then 'i' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[i, ]
              else if (i <= p+h) ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i-p]
              else { ## then 'i' refers to a SNP genotype (discrete)
                x <- as(Xsml[[sum(cumsum_sByChr < i-p-h)]][, i-p-h-cumsum_sByChr[sum(cumsum_sByChr < i-p-h)] ], "numeric")[, 1]+1
                if (any(x > 3, na.rm=TRUE))
                  warning(sprintf("i=%s has uncertain genotype calls which are treated here as missing", nam_i))
                x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
                nLevels[1] <- gLevels
              }
            }
            i <- 1L
            names(i) <- nam_i

            if (!is.na(nLevels[1]) || is.character(x) || is.factor(x) || is.logical(x)) {
              x <- factor(x)
              if (is.na(nLevels[1]))
                nLevels[1] <- nlevels(x)

              I <- 1L
            } else
              Y <- 1L

            missingMask[1] <- any(is.na(x))
            Xsub[, 1] <- as.numeric(x)

            nam_j <- j
            if (is.character(j)) {
              smt <- match(j, sNames)
              if (is.na(match(j, fNames)) && is.na(match(j, pNames)) && is.na(smt))
                stop(sprintf("j=%s does not form part of the variable names of the data\n", j))
              if (!is.na(match(j, fNames))) ## then 'j' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[j, ]
              else if (!is.na(match(j, pNames))) ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j]
              else { ## then 'j' refers to a SNP genotype (discrete)
                x <- as(Xsml[[sum(cumsum_sByChr < smt)]][, j], "numeric")+1
                if (any(x > 3, na.rm=TRUE))
                  warning(sprintf("j=%s has uncertain genotype calls which are treated here as missing", nam_j))
                x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
                nLevels[2] <- gLevels
              }
            } else {
              if (j <= p) ## then 'j' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[j, ]
              else if (j <= p+h) ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j-p]
              else { ## then 'j' refers to a SNP genotype (discrete)
                x <- as(Xsml[[sum(cumsum_sByChr < j-p-h)]][, j-p-h-cumsum_sByChr[sum(cumsum_sByChr < j-p-h)] ], "numeric")[, 1]+1
                if (any(x > 3, na.rm=TRUE))
                  warning(sprintf("j=%s has uncertain genotype calls which are treated here as missing", nam_j))
                x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call
                nLevels[2] <- gLevels
              }
            }
            j <- 2L
            names(j) <- nam_j

            if (!is.na(nLevels[2]) || is.character(x) || is.factor(x) || is.logical(x)) {
              x <- factor(x)
              if (is.na(nLevels[2]))
                nLevels[2] <- nlevels(x)
              I <- c(I, 2L)
            } else
              Y <- c(Y, 2L)

            missingMask[2] <- any(is.na(x))
            Xsub[, 2] <- as.numeric(x)

            nam_Q <- Q
            if (is.character(Q)) {
              Qe <- match(Q, fNames)
              if (any(!is.na(Qe))) {
                Xsub[, 3:(2+sum(!is.na(Qe)))] <- t(Biobase::exprs(X)[Qe[!is.na(Qe)], ])
                Y <- c(Y, 3:(2+sum(!is.na(Qe))))
              }
              if (any(is.na(Qe))) { ## then some variables in Q should refer to either phenotypic vars. or genotype calls
                Qp <- match(Q[is.na(Qe)], pNames)
                Qs <- c()
                if (any(is.na(Qp))) {
                  Qs <- match(Q[is.na(Qe)][is.na(Qp)], sNames)
                  if (any(is.na(Qs))) ## then some variables in Q should refer to genotype calls
                    stop(sprintf("Q={%s} do(es) not form part of the variable names of the data\n",
                                 paste(Q[is.na(Qe)][is.na(Qp)][is.na(Qs)], collapse=", ")))
                  Qp <- Qp[!is.na(Qp)]
                }
                for (k in seq(along=Qp)) {
                  x <- Biobase::pData(X)[, Qp[k]]
                  idx <- 2L+sum(!is.na(Qe))+k
                  missingMask[idx] <- any(is.na(x))

                  if (is.character(x) || is.factor(x) || is.logical(x)) {
                    x <- factor(x)
                    nLevels[idx] <- nlevels(x)
                    Xsub[, idx] <- as.numeric(x)
                    I <- c(I, idx)
                  } else {
                    Xsub[, idx] <- as.numeric(x)
                    Y <- c(Y, idx)
                  }
                }
                for (k in seq(along=Qs)) {
                  x <- as(Xsml[[sum(cumsum_sByChr < Qs[k])]][, Qs[k]-cumsum_sByChr[sum(cumsum_sByChr < Qs[k])] ], "numeric")[, 1]+1
                  if (any(x > 3, na.rm=TRUE))
                    warning(sprintf("Q=%s has uncertain genotype calls which are treated here as missing", Q[Qs[k]]))
                  x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call

                  idx <- 2L+sum(!is.na(Qe))+sum(!is.na(Qp))+k
                  missingMask[idx] <- any(is.na(x))
                  Xsub[, idx] <- x
                  I <- c(I, idx)
                  nLevels[idx] <- gLevels
                }
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            } else if (!is.null(Q)) { ## if argument Q was empty, it should remain empty
              if (any(Q > p+h+s))
                stop(sprintf("Q={%s} is/are larger than the number of expression profiles (%d), phenotypic variables (%d) and genotype calls (%d) together (%d+%d+%d=%d)\n",
                             paste(Q[Q > p+h+s], collapse=", "), p, h, s, p, h, s, p+h+s))

              if (any(Q <= p)) { ## Q indices smaller or equal than p correspond to expression profiles
                Xsub[, 3:(2+sum(Q <= p))] <- Biobase::exprs(X)[Q[Q <= p], ]
                Y <- c(Y, 3:(2+sum(Q <= p)))
              }
              Qp <- which(Q > p & Q <= p+h) ## Q indices larger than p and smaller than h correspond to phenotypic variables
              for (k in seq(along=Qp)) {
                x <- Biobase::pData(X)[, Q[Qp[k]]-p]
                idx <- 2L+sum(Q <= p)+k
                missingMask[idx] <- any(is.na(x))

                if (is.character(x) || is.factor(x)) {
                  x <- factor(x)
                  nLevels[idx] <- nlevels(x)
                  Xsub[, idx] <- as.numeric(x)
                  I <- c(I, idx)
                } else {
                  Xsub[, idx] <- as.numeric(x)
                  Y <- c(Y, idx)
                }
              }
              Qs <- which(Q > p+h) ## Q indices larger than p+h correspond to genotype calls
              for (k in seq(along=Qs)) {
                x <- as(Xsml[[sum(cumsum_sByChr < Q[Qs[k]]-p-h)]][, Q[Qs[k]]-p-h-cumsum_sByChr[sum(cumsum_sByChr < Q[Qs[k]]-p-h)] ], "numeric")[, 1]+1
                if (any(x > 3, na.rm=TRUE))
                  warning(sprintf("Q=%s has uncertain genotype calls which are treated here as missing", Q[Qs[k]]))
                x[x > 3] <- NA ## > 2 in the "numeric" coercion implies an uncertain call

                idx <- 2L+sum(Q <= p+h)+k
                missingMask[idx] <- any(is.na(x))
                Xsub[, idx] <- x
                I <- c(I, idx)
                nLevels[idx] <- gLevels
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            }

            rval <- NA

            if (is.null(I)) {
              S <- qpCov(Xsub)
              rval <- qpgraph:::.qpCItest(S, n, i, j, Q, R.code.only)
            } else {
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(Xsub)[I[nLevels[I]==1]]))

              missingData <- any(missingMask)
              ssd <- mapX2ssd <- NULL
              if (!missingData) {
                ssd <- qpCov(Xsub[, Y, drop=FALSE], corrected=FALSE)
                mapX2ssd <- match(colnames(Xsub), colnames(ssd))
                names(mapX2ssd) <- colnames(Xsub)
              }

              rval <- qpgraph:::.qpCItestHMGM(Xsub, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                              exact.test, R.code.only)
              if (is.nan(rval$statistic))
                warning(paste(sprintf("CI test unavailable for i=%s, j=%s and Q={",
                                      i, j, paste(Q, collapse=", ")),
                                      "}. Try a smaller Q or increase n if you can\n"))
            }

            class(rval) <- "htest" ## this is kind of redundant but otherwise
                                   ## the object returned by the C function does
                                   ## not print by default as an 'htest' object

            return(rval)
          })

## X comes as an ExpressionSet object
setMethod("qpCItest", signature(X="ExpressionSet"),
          function(X, i=1, j=2, Q=c(), exact.test=TRUE, R.code.only=FALSE) {
            p <- as.integer(nrow(X))
            h <- as.integer(ncol(pData(X)))
            n <- as.integer(ncol(X))
            fNames <- Biobase::featureNames(X)
            pNames <- colnames(Biobase::pData(X))
            Xsub <- matrix(0, nrow=n, ncol=length(c(i, j, Q)))
            colnames(Xsub) <- 1:length(c(i, j, Q))
            x <- Y <- I <- c()
            missingMask <- rep(FALSE, length(c(i, j, Q)))

            nam_i <- i
            if (is.character(i)) {
              if (is.na(match(i, fNames)) && is.na(match(i, pNames)))
                stop(sprintf("i=%s does not form part of the variable names of the data\n", i))
              if (!is.na(match(i, fNames))) ## then 'i' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[i, ]
              else ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i]
            } else {
              if (i <= p) ## then 'i' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[i, ]
              else ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, i-p]
            }
            i <- 1L
            names(i) <- nam_i

            if (is.character(x) || is.factor(x) || is.logical(x)) {
              Xsub[, 1] <- as.numeric(factor(x, levels=unique(x)))
              I <- 1L
            } else {
              Xsub[, 1] <- as.numeric(x)
              Y <- 1L
            }
            missingMask[1] <- any(is.na(x))

            nam_j <- j
            if (is.character(j)) {
              if (is.na(match(j, fNames)) && is.na(match(j, pNames)))
                stop(sprintf("j=%s does not form part of the variable names of the data\n", j))
              if (!is.na(match(j, fNames))) ## then 'j' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[j, ]
              else ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j]
            } else {
              if (j <= p) ## then 'j' refers to an expression profile (cont.)
                x <- Biobase::exprs(X)[j, ]
              else ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- Biobase::pData(X)[, j-p]
            }
            j <- 2L
            names(j) <- nam_j

            if (is.character(x) || is.factor(x) || is.logical(x)) {
              Xsub[, 2] <- as.numeric(factor(x, levels=unique(x)))
              I <- c(I, 2L)
            } else {
              Xsub[, 2] <- as.numeric(x)
              Y <- c(Y, 2L)
            }
            missingMask[2] <- any(is.na(x))

            nam_Q <- Q
            if (is.character(Q)) {
              Qe <- match(Q, fNames)
              if (any(!is.na(Qe))) {
                Xsub[, 3:(2+sum(!is.na(Qe)))] <- t(Biobase::exprs(X)[Qe[!is.na(Qe)], ])
                Y <- c(Y, 3:(2+sum(!is.na(Qe))))
              }
              if (any(is.na(Qe))) { ## then some variables in Q should refer to phenotypic vars.
                Qp <- match(Q[is.na(Qe)], pNames)
                if (any(is.na(Qp)))
                  stop(sprintf("Q={%s} do(es) not form part of the variable names of the data\n",
                       paste(Q[is.na(Qe)][is.na(Qp)], collapse=", ")))
                for (k in seq(along=Qp)) {
                  x <- Biobase::pData(X)[, Qp[k]]
                  idx <- 2L+sum(!is.na(Qe))+k

                  if (is.character(x) || is.factor(x)) {
                    Xsub[, idx] <- as.numeric(factor(x, levels=unique(x)))

                    I <- c(I, idx)
                  } else {
                    Xsub[, idx] <- as.numeric(x)
                    Y <- c(Y, idx)
                  }
                  missingMask[idx] <- any(is.na(x))
                }
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            } else if (!is.null(Q)) { ## if argument Q was empty, it should remain empty
              if (any(Q > p+h))
                stop(sprintf("Q={%s} is/are larger than the number of expression profiles (%d) and phenotypic variables (%d) together (%d+%d=%d)\n",
                             paste(Q[Q > p+h], collapse=", "), p, h, p, h, p+h))

              if (any(Q <= p)) { ## Q indices smaller or equal than p correspond to expression profiles
                Xsub[, 3:(2+sum(Q <= p))] <- Biobase::exprs(X)[Q[Q <= p], ]
                Y <- c(Y, 3:(2+sum(Q <= p)))
              }
              Qp <- which(Q > p) ## Q indices larger than p correspond to phenotypic variables
              for (k in seq(along=Qp)) {
                x <- Biobase::pData(X)[, Q[Qp[k]]-p]
                idx <- 2L+sum(Q <= p)+k

                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  Xsub[, idx] <- as.numeric(factor(x, levels=unique(x)))
                  I <- c(I, idx)
                } else {
                  Xsub[, idx] <- as.numeric(x)
                  Y <- c(Y, idx)
                }
                missingMask[idx] <- any(is.na(x))
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            }

            rval <- NA

            if (is.null(I)) {
              S <- qpCov(Xsub)
              rval <- qpgraph:::.qpCItest(S, n, i, j, Q, R.code.only)
            } else {
              missingData <- any(missingMask)
              ssd <- mapX2ssd <- NULL
              if (!missingData) {
                ssd <- qpCov(Xsub[, Y, drop=FALSE], corrected=FALSE)
                mapX2ssd <- match(colnames(Xsub), colnames(ssd))
                names(mapX2ssd) <- colnames(Xsub)
              }

              nLevels <- rep(NA_integer_, times=ncol(Xsub))
              nLevels[I] <- apply(Xsub[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(Xsub)[I[nLevels[I]==1]]))

              rval <- qpgraph:::.qpCItestHMGM(Xsub, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                              exact.test, R.code.only)
              if (is.nan(rval$statistic))
                warning(paste(sprintf("CI test unavailable for i=%s, j=%s and Q={",
                                      i, j, paste(Q, collapse=", ")),
                                      "}. Try a smaller Q or increase n if you can\n"))
            }

            class(rval) <- "htest" ## this is kind of redundant but otherwise
                                   ## the object returned by the C function does
                                   ## not print by default as a 'htest' object

            return(rval)
          })

## X comes as a data frame
setMethod("qpCItest", signature(X="data.frame"),
          function(X, i=1, j=2, Q=c(), I=NULL, long.dim.are.variables=TRUE,
                   exact.test=TRUE, R.code.only=FALSE) {

            X <- as.matrix(X)
            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)

            nam_i <- names(i)
            if (is.character(i)) {
              if (is.na(match(i, colnames(X))))
                stop(sprintf("i=%s does not form part of the variable names of the data\n",i))
              i <- match(i,colnames(X))
            }
            i <- as.integer(i)
            names(i) <- nam_i

            nam_j <- names(j)
            if (is.character(j)) {
              if (is.na(match(j, colnames(X))))
                stop(sprintf("j=%s does not form part of the variable names of the data\n",j))
              j <- match(j,colnames(X))
            }
            j <- as.integer(j)
            names(j) <- nam_j

            nam_Q <- names(Q)
            if (is.character(Q)) {
              if (any(is.na(match(Q, colnames(X)))))
                stop(sprintf("%s in Q does not form part of the variable names of the data\n",
                     Q[is.na(match(Q, colnames(X)))]))
              Q <- match(Q, colnames(X))
            }
            Q <- as.integer(Q)
            names(Q) <- nam_Q

            if (is.character(I)) {
              if (any(is.na(match(I, colnames(X)))))
                stop(sprintf("%s in I does not form part of the variable names of the data\n",
                     I[is.na(match(I, colnames(X)))]))
              I <- match(I, colnames(X))
            }

            rval <- NA

            if (is.null(I)) {
              S <- qpCov(X)
              n <- nrow(X)

              rval <- qpgraph:::.qpCItest(S, n, i, j, Q, R.code.only)
            } else {
              if (!is.character(I) && !is.numeric(I) && !is.integer(I))
                stop("I should be either variables names or indices\n")

              Y <- colnames(X)
              if (is.character(I))
                Y <- setdiff(colnames(X), I)
              else
                Y <- (1:ncol(X))[-I]

              if (is.character(Y)) {
                if (any(is.na(match(Y, colnames(X)))))
                  stop(sprintf("Some variables in Y do not form part of the variable names in X\n",i))
                Y <- match(Y, colnames(X))
              }

              nLevels <- rep(NA_integer_, times=ncol(X))
              nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
              if (any(nLevels[I] == 1))
                stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

              missingMask <- apply(X[, I, drop=FALSE], 2, function(x) any(is.na(x)))
              missingData <- any(missingMask)
              ssd <- mapX2ssd <- NULL
              if (!missingData) {
                ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
                mapX2ssd <- match(colnames(X), colnames(ssd))
                names(mapX2ssd) <- colnames(X)
              }

              rval <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                              exact.test, R.code.only)
              if (is.nan(rval$statistic))
                warning(paste(sprintf("CI test unavailable for i=%s, j=%s and Q={",
                                      colnames(X)[i], colnames(X)[j]),
                              paste(colnames(X)[Q], collapse=", "),
                                    "}. Try a smaller Q or increase n if you can\n"))
            }

            class(rval) <- "htest" ## this is kind of redundant but otherwise
                                   ## the object returned by the C function does
                                   ## not print by default as a 'htest' object

            return(rval)
          })

          
## X comes as a matrix
setMethod("qpCItest", signature(X="matrix"),
          function(X, i=1, j=2, Q=c(), I=NULL, n=NULL,
                   long.dim.are.variables=TRUE, exact.test=TRUE,
                   R.code.only=FALSE) {

            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X)))
              colnames(X) <- 1:ncol(X)

            nam_i <- i
            if (is.character(i)) {
              if (is.na(match(i, colnames(X))))
                stop(sprintf("i=%s does not form part of the variable names of the data\n",i))
              i <- match(i,colnames(X))
            }
            i <- as.integer(i)
            names(i) <- nam_i

            nam_j <- j
            if (is.character(j)) {
              if (is.na(match(j, colnames(X))))
                stop(sprintf("j=%s does not form part of the variable names of the data\n",j))
              j <- match(j,colnames(X))
            }
            j <- as.integer(j)
            names(j) <- nam_j

            nam_Q <- Q
            if (is.character(Q)) {
              if (any(is.na(match(Q, colnames(X)))))
                stop(sprintf("%s in Q does not form part of the variable names of the data\n",
                     Q[is.na(match(Q, colnames(X)))]))
              Q <- match(Q, colnames(X))
            }
            Q <- as.integer(Q)
            names(Q) <- nam_Q

            if (is.character(I)) {
              if (any(is.na(match(I, colnames(X)))))
                stop(sprintf("%s in I does not form part of the variable names of the data\n",
                     I[is.na(match(I, colnames(X)))]))
              I <- match(I, colnames(X))
            }

            rval <- NA

            # if the matrix is squared let's assume then that it is the sample
            # covariance matrix and that the sample size is the next parameter
            if (nrow(X) != ncol(X)) {
              if (!is.null(n))
                stop("If X is not a sample covariance matrix then N should not be set\n")

              if (is.null(I)) {
                S <- qpCov(X)
                n <- nrow(X)

                rval <- qpgraph:::.qpCItest(S, n, i, j, Q, R.code.only)
              } else {
                if (!is.character(I) && !is.numeric(I) && !is.integer(I))
                  stop("I should be either variables names or indices\n")

                Y <- colnames(X)
                if (is.character(I))
                  Y <- setdiff(colnames(X), I)
                else
                  Y <- (1:ncol(X))[-I]

                if (is.character(Y)) {
                  if (any(is.na(match(Y, colnames(X)))))
                    stop(sprintf("Some variables in Y do not form part of the variable names in X\n",i))
                  Y <- match(Y, colnames(X))
                }

                nLevels <- rep(NA_integer_, times=ncol(X))
                nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
                if (any(nLevels[I] == 1))
                  stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

                missingMask <- apply(X[, I, drop=FALSE], 2, function(x) any(is.na(x)))
                missingData <- any(missingMask)
                ssd <- mapX2ssd <- NULL
                if (!missingData) {
                  ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
                  mapX2ssd <- match(colnames(X), colnames(ssd))
                  names(mapX2ssd) <- colnames(X)
                }

                rval <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                                exact.test, R.code.only)
                if (is.nan(rval$statistic))
                  warning(paste(sprintf("CI test unavailable for i=%s, j=%s and Q={",
                                        colnames(X)[i], colnames(X)[j]),
                                paste(colnames(X)[Q], collapse=", "),
                                      "}. Try a smaller Q or increase n if you can\n"))
              }

            } else {
              if (!is.null(I))
                stop("If X is a sample covariance matrix then I should not be set\n")

              if (is.null(rownames(X)))
                rownames(X) <- colnames(X)

              S <- X
              rval <- qpgraph:::.qpCItest(S, n, i, j, Q, R.code.only)
            }

            class(rval) <- "htest" ## this is kind of redundant but otherwise
                                   ## the object returned by the C function does
                                   ## not print by default as a 'htest' object

            return(rval)
          })

.qpCItest <- function(S, n, i=1L, j=2L, Q=c(), R.code.only=FALSE) {

  p <- (d <- dim(S))[1]
  if (p != d[2] || !isSymmetric(S))
    stop("S is not squared and symmetric. Is it really a covariance matrix?n")

  if (!is.integer(i) || !is.integer(j) || (!is.null(Q) && !is.integer(Q)))
    stop("i, j and Q should contain only integer values when calling .qpCItest()")

  if (!R.code.only) {
    return(qpgraph:::.qpFastCItestStd(S, n, i, j, Q));
  }

  q       <- length(Q)
  Mmar    <- S[c(i, j, Q), c(i, j, Q)]
  S11     <- Mmar[1,1]
  S12     <- Mmar[1,-1]
  S21     <- Mmar[-1,1]
  S22     <- Mmar[-1,-1]
  S22inv  <- solve(S22)
  betahat <- as.numeric(S12 %*% S22inv[,1])
  df      <- c(df = n - q - 2)
  sigma   <- sqrt((S11 - S12 %*% S22inv %*% S21) * (n - 1) / df)
  se      <- as.numeric(sigma * sqrt(S22inv[1,1] / (n - 1)))
  t.value <- c(t = betahat / se)
  p.value <- c(two.sided = 2 * (1 - pt(abs(t.value), df)))
  est     <- c(beta = betahat)
  n.value <- c("partial regression coefficient" = 0)
  method  <- "Conditional independence test for continuous data using a t test for zero partial regression coefficient"

  RVAL <- list(statistic=t.value,
               parameter=df,
               p.value=p.value,
               estimate=est,
               null.value=n.value,
               alternative="two.sided",
               method=method,
               data.name=sprintf("%s and %s given {%s}", names(i), names(j), paste(names(Q), collapse=", ")))
  class(RVAL) <- "htest"

  RVAL
}

.ssdStats <- function(X, I, Y) {

  n <- n_co <- dim(X)[1]

  if (length(I) == 0)
    return((n-1) * cov(X[, Y, drop=FALSE]))

  missingMask <- apply(X[, I, drop=FALSE], 1, function(x) any(is.na(x)))

  xtab <- tapply(1:n, as.data.frame(X[, I, drop=FALSE]))
  xtab[missingMask] <- -1 ## label missing observations
  xtab <- split(as.data.frame(X[, Y, drop=FALSE]), xtab)
  xtab <- xtab[as.integer(names(xtab)) > 0] ## remove missing observations
  xtab <- xtab[which(sapply(lapply(xtab, dim), "[", 1) > 0)]
  ni <- sapply(lapply(xtab, dim), "[", 1)
  ssd <- Reduce("+",
                lapply(as.list(1:length(xtab)),
                       function(i, x) qpCov(as.matrix(x[[i]]), corrected=FALSE), xtab))
                       ##function(i, x, ni, n) (ni[i]-1)*cov(x[[i]]), xtab, ni, n))
  ssd
}

.qpCItestHMGM <- function(X, I, nLevels, Y, ssdMat, mapX2ssdMat, i, j, Q,
                          exact.test=TRUE, R.code.only=FALSE ) {
  if (!is.null(ssdMat)) {
    p <- (d <- dim(ssdMat))[1]
    if (p != d[2] || !isSymmetric(ssdMat))
      stop("ssdMat is not squared and symmetric. Is it really an ssd matrix?\n")
  }

  ## not possible when considering restrict.Q != NULL
  ## if (p != length(Y))
  ##  stop("ssdMat is not the ssd matrix of the variables in Y\n")

  ## not necessary if we assume .processParameters() is called before
  ## if (!is.integer(i) || !is.integer(j) || (!is.null(Q) && !is.integer(Q)) || (!is.null(I) && !is.integer(I)) || !is.integer(Y))
  ##   stop("i, j, Y, I and Q should contain only integer values when calling .qpCItestHMGM()")

  if (all(!is.na(match(c(i,j), I))))
    stop("i and j cannot be both discrete at the moment")

  ## for the sake of speed assume the map is correct
  ## if (is.null(rownames(ssdMat)) || is.null(colnames(ssdMat)) ||
  ##     any(is.na(match(colnames(ssdMat), colnames(X)))))
  ##   stop("ssdMat should have row and column names that match variables in X\n")

  if (!is.na(match(j, I))) { ## if any of (i,j) is discrete, it should be i
    tmp <- i
    i <- j
    j <- tmp
    tmp <- nLevels[i]
    nLevels[i] <- nLevels[j]
    nLevels[j] <- tmp
  }

  if (!R.code.only) {
    return(qpgraph:::.qpFastCItestHMGM(X, I, nLevels, Y, ssdMat, mapX2ssdMat,
                                       i, j, Q, exact.test))
  }

  I <- intersect(I, c(i, Q))
  Y <- intersect(Y, c(i, j, Q))

  ssd <- ssd_i <- ssd_j <- ssd_ij <- diag(2) 
  n <- nrow(X)

  if (length(I) == 0) { ## dspMatrix -> matrix because det() still doesn't work with Matrix classes
    ## when we'll handle missing continuous data too we'll have to consider calculating ssd
    ssd <- as.matrix(ssdMat[mapX2ssdMat[Y], mapX2ssdMat[Y], drop=FALSE])
    ## ssd_i = ssd_Gamma when i is discrete or ssd_{Gamma\i} when i is continuous
    ssd_i <- as.matrix(ssdMat[mapX2ssdMat[setdiff(Y, i)], mapX2ssdMat[setdiff(Y, i)], drop=FALSE])
    if (length(setdiff(Y, j)) > 0) {
      ssd_j <- as.matrix(ssdMat[mapX2ssdMat[setdiff(Y, j)], mapX2ssdMat[setdiff(Y, j)], drop=FALSE])
      if (length(setdiff(Y, c(i, j))) > 0)
        ssd_ij <- as.matrix(ssdMat[mapX2ssdMat[setdiff(Y, c(i, j))],
                                   mapX2ssdMat[setdiff(Y, c(i, j))], drop=FALSE])
    }
  } else {
    ## by now missing data w/ complete.obs
    missingMask <- apply(X[, I, drop=FALSE], 1, function(x) any(is.na(x)))
    missingData <- any(missingMask)
    n_co <- n - sum(missingMask)
    ssd <- qpgraph:::.ssdStats(X[!missingMask, ], I, Y)
    if (length(setdiff(I, i)) == 0 && !missingData && !is.null(ssdMat)) ## dspMatrix -> matrix because det() still doesn't work with Matrix classes
      ssd_i <- as.matrix(ssdMat[mapX2ssdMat[setdiff(Y, i)], mapX2ssdMat[setdiff(Y, i)], drop=FALSE])
    else ## ssd_i = ssd_Gamma when i is discrete or ssd_{Gamma\i} when i is continuous
      ssd_i <- qpgraph:::.ssdStats(X[!missingMask, ], setdiff(I, i), setdiff(Y, i))

    if (length(setdiff(Y, j)) > 0) {
      ssd_j <- qpgraph:::.ssdStats(X[!missingMask, ], I, setdiff(Y, j))
      if (length(setdiff(Y, c(i,j))) > 0) {
        if (length(setdiff(I, i)) == 0 && !missingData && !is.null(ssdMat)) ## dspMatrix -> matrix because det() still doesn't work with Matrix classes
          ssd_ij <- as.matrix(ssdMat[mapX2ssdMat[setdiff(Y, c(i, j))],
                                     mapX2ssdMat[setdiff(Y, c(i, j))], drop=FALSE])
        else
          ssd_ij <- qpgraph:::.ssdStats(X[!missingMask, ], setdiff(I, i), setdiff(Y, c(i, j)))
      }
    }
  }
 
  ssd <- determinant(ssd, logarithm=TRUE)
  ssd_ij <- determinant(ssd_ij, logarithm=TRUE)
  ssd_j <- determinant(ssd_j, logarithm=TRUE)
  ssd_i <- determinant(ssd_i, logarithm=TRUE)
  final_sign <- ssd$sign * ssd_ij$sign * ssd_j$sign *ssd_i$sign

  ## lr <- -n * log((det(ssd) * det(ssd_ij)) / 
  ##               (det(ssd_j) * det(ssd_i)))

  lr <- NaN
  p.value <- df <- a <- b <- NA
  stat <- param <- n.value <- method <- alt <- NA

  zero_boundary <- log(.Machine$double.eps)
  if (ssd$modulus > zero_boundary && ssd_ij$modulus > zero_boundary &&
      ssd_j$modulus > zero_boundary && ssd_i$modulus > zero_boundary &&
      final_sign == 1) {
    mixedEdge <- sum(!is.na(match(c(i, j), I))) > 0
    nGamma <- length(Y)
    Delta <- I
    DeltaStar <- setdiff(I, c(i, j))
    if (exact.test) {
      lr <- exp(ssd$modulus[1]+ssd_ij$modulus[1]-ssd_j$modulus[1]-ssd_i$modulus[1])
      a <- (n_co-nGamma-prod(nLevels[Delta])+1)/2 ## by now missing data w/ complete.obs
      b <- ifelse(mixedEdge,
                  prod(nLevels[DeltaStar])*(nLevels[intersect(Delta, c(i,j))]-1)/2,
                  0.5)
      if (a > 0 && b > 0)
        p.value <- c(less = pbeta(q=lr, shape1=a, shape2=b, lower.tail=TRUE))
      else {
        p.value <- a <- b <- NA
        names(p.value) <- "less"
      }
      stat <- lr
      names(stat) <- "Lambda"
      param <- c(a=a, b=b)
      n.value <- c("Lambda" = 1)
      method  <- "Conditional independence test for homogeneous mixed data using an exact likelihood ratio test"
      alt <- "less"
    } else {
      lr <- -n_co * (ssd$modulus[1]+ssd_ij$modulus[1]-ssd_j$modulus[1]-ssd_i$modulus[1])
      df <- 1
      if (mixedEdge)
        df <- prod(nLevels[DeltaStar])*(nLevels[intersect(Delta, c(i, j))]-1)
      stat <- lr
      names(stat) <- "-n log Lambda"
      param <- c(df=df)
      p.value <- c(greater = 1 - pchisq(lr, df=df, lower.tail=TRUE))
      n.value <- c("-n log Lambda" = 0)
      method  <- "Conditional independence test for homogeneous mixed data using an asymptotic likelihood ratio test"
      alt <- "greater"
    }
  }

  RVAL <- list(statistic=stat,
               parameter=param,
               p.value=p.value,
               estimate=NULL,
               null.value=n.value,
               alternative=alt,
               method=method,
               data.name=sprintf("%s and %s given {%s}", names(i), names(j), paste(names(Q), collapse=", ")))
  class(RVAL) <- "htest"

  RVAL
}


.qpFastCItestStd <- function(S, n, i, j, Q) {
  return(.Call("qp_fast_ci_test_std", S@x, nrow(S), as.integer(n), i, j, Q))
}

.qpFastCItestOpt <- function(S, n, i, j, Q) {
  return(.Call("qp_fast_ci_test_opt", S@x, nrow(S), as.integer(n), i, j, Q))
}

.qpFastCItestHMGM <- function(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q, exact.test) {
  x <- NULL
  if (!is.null(ssd))
    x <- ssd@x
  return(.Call("qp_fast_ci_test_hmgm", X, I, nLevels, Y, x, as.integer(mapX2ssd),
               i, j, Q, as.integer(exact.test)))
}
