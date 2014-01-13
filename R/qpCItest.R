## qpgraph package - this R code implements functions to learn qp-graphs from
##                   data, to estimate partial correlations, simulate undirected Gaussian
##                   graphical models and deal with microarray and genetic data in order
##                   to build network models of molecular regulation
##
## Copyright (c) 2008-2013 R. Castelo and A. Roverato, with contributions from Inma Tur.
## This package is open source and free software; you can redistribute it and/or
## modify it under the terms of the Artistic License 2.0
## as published at http://www.r-project.org/Licenses/Artistic-2.0
##
## Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT
## HOLDER AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED
## WARRANTIES.  THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
## PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT
## PERMITTED BY YOUR LOCAL LAW.  UNLESS REQUIRED BY LAW, NO COPYRIGHT
## HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT,
## INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE
## OF THE PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## See the Artistic License 2.0 for more details.


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

## ADD gLevels to either directly denote the number of genotype levels or trigger their
## automatic calculation

## X comes as a GGBase::smlSet object
setMethod("qpCItest", signature(X="smlSet"),
          function(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

            p <- as.integer(nrow(X))
            h <- as.integer(ncol(Biobase::pData(X)))
            Xsml <- GGBase::smList(X)
            sByChr <- sapply(Xsml, ncol)
            gLevels <- sum(unique(as.vector(as(Xsml[[1]][, 1:min(sByChr[1], 1000)], "matrix"))) > 0, na.rm=TRUE)
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
            rownames(Xsub) <- 1:nrow(Xsub)

            if (is.null(I)) {
              S <- qpCov(Xsub)
              rval <- qpgraph:::.qpCItest(S, i, j, Q, R.code.only)
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
                                              exact.test, use, tol, R.code.only)
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

## X comes as a Biobase::ExpressionSet object
setMethod("qpCItest", signature(X="ExpressionSet"),
          function(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

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
            rownames(Xsub) <- 1:nrow(Xsub)

            if (is.null(I)) {
              S <- qpCov(Xsub)
              rval <- qpgraph:::.qpCItest(S, i, j, Q, R.code.only)
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
                                              exact.test, use, tol, R.code.only)
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

## X comes as a qtl::cross object
setMethod("qpCItest", signature(X="cross"),
          function(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                   tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

            p <- as.integer(qtl::nphe(X))
            sByChr <- qtl::nmar(X)
            gLevels <- sum(unique(as.vector(as(X$geno[[1]]$data[, 1:min(sByChr[1], 1000)], "matrix"))) > 0, na.rm=TRUE)
            cumsum_sByChr <- c(0, cumsum(sByChr))
            s <- sum(sByChr)
            n <- as.integer(qtl::nind(X))
            pNames <- colnames(X$pheno)
            sNamesByChr <- lapply(X$geno, function(x) colnames(x$data))
            sNames <- unlist(sNamesByChr, use.names=FALSE)
            qtlNames <- c()
            nqtl <- 0 
            if ("qtlgeno" %in% names(X)) {
              qtlNames <- colnames(X$qtlgeno)
              nqtl <- length(qtlNames)
            }
            Xsub <- matrix(0, nrow=n, ncol=length(c(i, j, Q)))
            colnames(Xsub) <- 1:length(c(i, j, Q))
            nLevels <- rep(NA, length(c(i, j, Q)))
            missingMask <- rep(FALSE, length(c(i, j, Q)))
            x <- Y <- I <- c()

            nam_i <- i
            if (is.character(i)) {
              smt <- match(i, sNames)
              if (is.na(match(i, qtlNames)) && is.na(match(i, pNames)) && is.na(smt))
                stop(sprintf("i=%s does not form part of the variable names of the data\n", i))
              if (!is.na(match(i, pNames))) ## then 'i' refers to a phenotypic variable (cont. or discrete)
                x <- X$pheno[, i]
              else if (!is.na(match(i, qtlNames))) { ## then 'i' refers to a QTL genotype (discrete)
                x <- as(X$qtlgeno[, i], "numeric") ## assume genotypes come as 1, 2, ...
                nLevels[1] <- gLevels
              } else { ## then 'i' refers to a marker genotype (discrete)
                x <- as(X$geno[[sum(cumsum_sByChr < smt)]]$data[, i], "numeric") ## assume genotypes come as 1, 2, ...
                nLevels[1] <- gLevels
              }
            } else {
              if (i <= p) ## then 'i' refers to a phenotype
                x <- X$pheno[, i]
              else if (i <= p+s) ## then 'i' refers to a marker genotype (discrete)
                x <- as(X$geno[[sum(cumsum_sByChr < i-p)]]$data[, i-p], "numeric") ## assume genotypes come as 1, 2, ...
              else { ## then 'i' refers to a QTL (discrete)
                x <- as(X$qtlgeno[, i-p-s], "numeric") ## assume genotypes come as 1, 2, ...
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
              if (is.na(match(j, qtlNames)) && is.na(match(j, pNames)) && is.na(smt))
                stop(sprintf("j=%s does not form part of the variable names of the data\n", j))
              if (!is.na(match(j, pNames))) ## then 'j' refers to a phenotypic variable (cont. or discrete)
                x <- X$pheno[, j]
              else if (!is.na(match(j, qtlNames))) { ## then 'i' refers to a QTL genotype (discrete)
                x <- as(X$qtlgeno[, j], "numeric") ## assume genotypes come as 1, 2, ...
                nLevels[2] <- gLevels
              } else { ## then 'j' refers to a marker genotype (discrete)
                x <- as(X$geno[[sum(cumsum_sByChr < smt)]]$data[, j], "numeric") ## assume genotypes come as 1, 2, ...
                nLevels[2] <- gLevels
              }
            } else {
              if (j <= p) ## then 'j' refers to a phenotype
                x <- X$pheno[, j]
              else if (j <= p+s) ## then 'j' refers to a marker genotype (discrete)
                x <- as(X$geno[[sum(cumsum_sByChr < j-p)]]$data[, j-p], "numeric") ## assume genotypes come as 1, 2, ...
              else { ## then 'j' refers to a QTL (discrete)
                x <- as(X$qtlgeno[, j-p-s], "numeric") ## assume genotypes come as 1, 2, ... 
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
              Qp <- match(Q, pNames)
              Qs <- Qqtl <- c()
              if (any(is.na(Qp))) { ## then some variables in Q should refer to marker genotypes
                Qs <- match(Q[is.na(Qp)], sNames)
                if (any(is.na(Qs))) { ## then some variables in Q should refer to QTL genotypes
                  Qqtl <- match(Q[is.na(Qp)][is.na(Qs)], qtlNames)
                  if (any(is.na(Qqtl)))
                    stop(sprintf("Q={%s} do(es) not form part of the variable names of the data\n",
                                 paste(Q[is.na(Qp)][is.na(Qs)][is.na(Qqtl)], collapse=", ")))
                  Qs <- Qs[!is.na(Qs)]
                }
                Qp <- Qp[!is.na(Qp)]
              }
              for (k in seq(along=Qp)) {
                x <- X$pheno[, Qp[k]]
                idx <- 2L+k
                missingMask[idx] <- any(is.na(x))

                if (is.character(x) || is.factor(x) || is.logical(x)) {
                  x <- factor(x)
                  nLevels[idx] <- nlevels(x)
                  Xsub[, idx] <- as.numeric(x) ## transforming factor variable values into 1, 2, ...
                  I <- c(I, idx)
                } else {
                  Xsub[, idx] <- as.numeric(x)
                  Y <- c(Y, idx)
                }
              }
              for (k in seq(along=Qs)) {
                x <- as(X$geno[[sum(cumsum_sByChr < Qs[k])]]$data[, Qs[k]-cumsum_sByChr[sum(cumsum_sByChr < Qs[k])] ], "numeric") ## assume genotypes come as 1, 2, ...
                idx <- 2L+sum(!is.na(Qp))+k
                missingMask[idx] <- any(is.na(x))
                Xsub[, idx] <- x
                I <- c(I, idx)
                nLevels[idx] <- gLevels
              }
              for (k in seq(along=Qqtl)) {
                x <- as(X$qtlgeno[, Qqtl[k]], "numeric") ## assume genotypes come as 1, 2, ...
                idx <- 2L+sum(!is.na(Qp))+sum(!is.na(Qs))+k
                missingMask[idx] <- any(is.na(x))
                Xsub[, idx] <- x
                I <- c(I, idx)
                nLevels[idx] <- gLevels
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            } else if (!is.null(Q)) { ## if argument Q was empty, it should remain empty
              if (any(Q > p+s+nqtl))
                stop(sprintf("Q={%s} is/are larger than the number of phenotypic variables (%d), marker genotypes (%d) and QTL genotypes together (%d+%d+%d=%d)\n",
                             paste(Q[Q > p+s+nqtl], collapse=", "), p, s, nqtl, p, s, nqtl, p+s+nqtl))

              Qp <- which(Q <= p) ## Q smaller than p correspond to phenotypic variables
              for (k in seq(along=Qp)) {
                x <- X$pheno[, Q[Qp[k]]]
                idx <- 2L+k
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
              Qs <- which(Q > p && Q <= p+s) ## Q indices larger than p and smaller than p+s correspond to marker genotypes
              for (k in seq(along=Qs)) {
                x <- as(X$geno[[sum(cumsum_sByChr < Q[Qs[k]]-p)]]$data[, Q[Qs[k]]-p-cumsum_sByChr[sum(cumsum_sByChr < Q[Qs[k]]-p)] ], "numeric") ## assume genotypes come as 1, 2, ...

                idx <- 2L+sum(Q <= p)+k
                missingMask[idx] <- any(is.na(x))
                Xsub[, idx] <- x
                I <- c(I, idx)
                nLevels[idx] <- gLevels
              }
              Qqtl <- which(Q > p+s) ## Q indices larger than p+s correspond to QTL genotypes
              for (k in seq(along=Qqtl)) {
                x <- as(X$qtlgeno[, Q[Qqtl[k]]-p-s], "numeric") ## assume gentypes come as 1, 2, ...
                idx <- 2L+sum(Q <= p+s)+k
                missingMask[idx] <- any(is.na(x))
                Xsub[, idx] <- x
                I <- c(I, idx)
                nLevels[idx] <- gLevels
              }
              Q <- 3:length(c(i, j, Q))
              names(Q) <- nam_Q
            }

            rval <- NA
            rownames(Xsub) <- 1:nrow(Xsub)

            if (is.null(I)) {
              S <- qpCov(Xsub)
              rval <- qpgraph:::.qpCItest(S, i, j, Q, R.code.only)
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
                                              exact.test, use, tol, R.code.only)
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



## X comes as a data frame
setMethod("qpCItest", signature(X="data.frame"),
          function(X, i=1, j=2, Q=c(), I=NULL, long.dim.are.variables=TRUE,
                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

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
            rownames(Xsub) <- 1:nrow(Xsub)

            if (is.null(I)) {
              S <- qpCov(X)
              n <- nrow(X)

              rval <- qpgraph:::.qpCItest(S, i, j, Q, R.code.only)
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
                                              exact.test, use, tol, R.code.only)
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
          function(X, i=1, j=2, Q=c(), I=NULL,
                   long.dim.are.variables=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE) {

            use <- match.arg(use)

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
            rownames(X) <- 1:nrow(X)

            if (is.null(I)) {
              if (use == "em")
                stop("EM not implemented yet for missing values in continuous variables. Please set use=\"complete.obs\"\n")

              S <- qpCov(X)
              n <- nrow(X)

              rval <- qpgraph:::.qpCItest(S, i, j, Q, R.code.only)
            } else {
              if (!is.character(I) && !is.numeric(I) && !is.integer(I))
                stop("argument I should contain either variables names or indices\n")

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

              missingMask <- apply(X, 1, function(x) any(is.na(x)))
              missingData <- any(missingMask)
              ssd <- mapX2ssd <- NULL
              if (!missingData) {
                ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
                ## mapX2ssd <- match(colnames(X), colnames(ssd))
                ## names(mapX2ssd) <- colnames(X)
                mapX2ssd <- rep(NA, ncol(X))
                mapX2ssd[Y] <- 1:length(Y)
              }

              rval <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                              exact.test, use, tol, R.code.only)
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

## X comes as an SsdMatrix (i.e., a covariance matrix calculated with qpCov())
setMethod("qpCItest", signature(X="SsdMatrix"),
          function(X, i=1, j=2, Q=c(), R.code.only=FALSE) {

            use <- match.arg(use)

            if (!is.double(X))
              stop("X should be double-precision real numbers\n")

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

            if (is.null(rownames(X)))
              rownames(X) <- colnames(X)

            rval <- qpgraph:::.qpCItest(X, i, j, Q, R.code.only)

            class(rval) <- "htest" ## this is kind of redundant but otherwise
                                   ## the object returned by the C function does
                                   ## not print by default as a 'htest' object

            return(rval)
          })

.qpCItest <- function(S, i=1L, j=2L, Q=c(), R.code.only=FALSE) {

  if (class(S) != "SsdMatrix")
    stop("internal function qpgraph:::.qpCItest() expects an 'SsdMatrix' object as first argument\n");

  p <- (d <- dim(S))[1]
  if (p != d[2] || !isSymmetric(S))
    stop("S is not squared and symmetric. Is it really a covariance matrix?n")
  n <- S@n

  if (!is.integer(i) || !is.integer(j) || (!is.null(Q) && !is.integer(Q)))
    stop("i, j and Q should contain only integer values when calling .qpCItest()")

  if (!R.code.only) {
    return(qpgraph:::.qpFastCItestStd(S, i, j, Q));
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

##
## calculate the ssd matrix using complete observations only
##
.ssdStatsCompleteObs <- function(X, I, Y, missingMask) {

  if (length(I) == 0) {
    ssd <- qpCov(X[!missingMask, Y, drop=FALSE], corrected=FALSE)
    return(ssd)
  }

  n <- nrow(X)
  n.co <- n - sum(missingMask)

  xtab <- tapply(1:n, as.data.frame(X[, I, drop=FALSE]))
  xtab[missingMask] <- -1 ## label missing observations
  xtab <- split(as.data.frame(X[, Y, drop=FALSE]), xtab)
  xtab <- xtab[as.integer(names(xtab)) > 0] ## remove missing observations
  xtab <- xtab[which(sapply(lapply(xtab, dim), "[", 1) > 0)]
  ## ni <- sapply(lapply(xtab, dim), "[", 1)
  ssd <- Reduce("+",
                lapply(as.list(1:length(xtab)),
                       function(i, x) qpCov(as.matrix(x[[i]]), corrected=FALSE)@ssd, xtab))
                        ##function(i, x, ni, n) (ni[i]-1)*cov(x[[i]]), xtab, ni, n))
  new("SsdMatrix", ssd=as(ssd, "dspMatrix"), n=n.co)
}

##
## the following functions calculate the ssd matrix using the EM algorithm
##

## k(i) = y^T\Sigma^{-1}\mu(i) - 1/2 * [y^T\Sigma^{-1} y + \mu(i)^T \Sigma^{-1}\mu(i)] + log p(i)
Ki <- function(x, Ys, i, mapX2Y, Sigma, mu, p) {

  aux <- p[i]
  if (length(Ys) > 0) {
    Ys_Sigma <- mapX2Y[Ys]
    aux <- exp( t(x[ ,Ys]) %*% solve(Sigma[Ys_Sigma, Ys_Sigma]) %*% mu[i, Ys_Sigma] -
                (t(x[ ,Ys]) %*% solve(Sigma[Ys_Sigma, Ys_Sigma]) %*% x[ ,Ys] +
                 t(mu[i, Ys_Sigma]) %*% solve(Sigma[Ys_Sigma, Ys_Sigma]) %*% mu[i, Ys_Sigma])/2 + log(p[i]) )
  }

  aux
}

## calculate the probability of I=i for each observation with missing values
## pr(I=i' | (i_{obs}, y)^{(\nu)}) = exp k(i') / \sum_{s\in{\cal S}} exp k(s)
prob_i <- function(idxMissingObs, X, I, Is, Ys, k, level_i, levels_I, mapX2I, mapX2Y, p, mu, Sigma) {

  p_i <- vector(length=length(idxMissingObs))
  ## names(p_i) <- idxMissingObs

  for (i in 1:length(idxMissingObs)) {
    x <- X[idxMissingObs[i], ,drop=FALSE]
    ## I_obs <- which(!is.na(x[, I])) ## 2/4/13 I_obs should be on the X-scale and not on the I-scale
    I_obs <- intersect(which(!is.na(x)), I)

    ## if (length(I_obs) > 0 && !all(x[, intersect(Is, I_obs)] == level_i[intersect(Is, I_obs)])) {
    if (length(I_obs) > 0 && !all(x[, intersect(Is, I_obs)] == level_i[!is.na(match(Is, I_obs))])) {
      p_i[i] <- 0
    } else {
      ifelse(length(I_obs)==0,
             index_S <- 1:nrow(levels_I),
             index_S <- which(apply(levels_I[, mapX2I[I_obs], drop=FALSE], 1, function(l) {all(l == x[, I_obs])})))
      index_Si <- index_S[which(sapply(index_S, function(i) {all(levels_I[i, mapX2I[Is]] == level_i)}))]
      if (all(intersect(index_S, index_Si)==index_S)) {
        p_i[i] <- 1
      } else {
        K_den <- sum(sapply(index_S, function(i) {Ki(x, Ys, i, mapX2Y, Sigma, mu, p)}))
        ifelse (K_den == 0,
                p_i[i] <- 0,
                p_i[i] <- sum(sapply(index_Si, function(i) {Ki(x, Ys, i, mapX2Y, Sigma, mu, p)/K_den})))
      }
    }
  }

  p_i
}

## complete sufficient statistics: calculate sufficient statistics (En, Es and Ess) of those
## observations that do not have missing values
stat_com <- function(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs, Is, Ys, levels_Is) {
  Es_com <- Ess_com <- n_com <- c()

  if (length(idxCompleteObs) > 0) {
    xi_level <- lapply(1:nrow(levels_Is), function(k) names(which(apply(X[idxCompleteObs, Is, drop=FALSE], 1, function(x){all(x == levels_Is[k, ])}))))
    n_com <- sapply(xi_level, length)
  }

  if (length(Ys) > 0) {
    Es_com <- matrix(0, nrow=nrow(levels_Is), ncol=length(Ys))##, dimnames=list(1:nrow(levels_Is), Ys))
    Ess_com <- matrix(0, nrow=length(Ys), ncol=length(Ys))##, dimnames=list(Ys, Ys))
    if (length(idxCompleteObs) > 0) {
      aux <- which(n_com > 0)
      Es_com[aux, ] <- matrix(t(sapply(aux, function(i) colSums(X[xi_level[[i]], Ys, drop=FALSE]))), ncol=length(Ys))##, dimnames=list(aux, Ys))
      Ess_com <- Reduce("+" ,lapply(xi_level, function(i) t(as.matrix(X[i, Ys, drop=FALSE]))%*%as.matrix(X[i, Ys,drop=FALSE])))
    }
  }

  list(idxMissingObs=idxMissingObs, mapAllObs2MissingObs=mapAllObs2MissingObs,
       n_com=n_com, Es_com=Es_com, Ess_com=Ess_com)
}

## sufficient statistics from missing data: perform the E step for the observations with missing values.
stat_mis <- function(X, Is, Ys, levels_Is, stat, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma) {

  m <- vector(length=nrow(levels_I))
  h <- bar_y <- matrix(0, nrow=nrow(levels_I), ncol=length(Y))##, dimnames=list(1:nrow(levels_I), Y))
  K <- matrix(0, nrow=length(Y), ncol=length(Y))##, dimnames=list(Y, Y))
  n <- nrow(X)

  Ys_Sigma <- mapX2Y[Ys]

  ### Is != emptyset
  if (length(Is) > 0) {

    ### Pr(I=i | x_obs)
    p_i <- sapply(1:nrow(levels_Is),
                  function(k) prob_i(stat$idxMissingObs, X, I, Is, Ys, k,
                                     levels_Is[k, ], levels_I, mapX2I, mapX2Y, p, mu, Sigma))
    if (!is.matrix(p_i)) {
      p_i <- t(as.matrix(p_i))
      ## dimnames(p_i) <- list(stat$idxMissingObs, 1:nrow(levels_Is))
    }

    index <- lapply(1:nrow(levels_Is), function(k) which(apply(levels_I[, mapX2I[Is], drop=FALSE], 1, function(l){all(l == levels_Is[k, ])})))

    ### En
    E_n <- stat$n_com + colSums(p_i)

    ### m     # afegit
    for (k in 1:nrow(levels_Is)) {
      j <- index[[k]]
      m[j] <- rep(E_n[k], length(j))
    }

    ### Is != emptyset & Ys != emptyset
    if (length(Ys) > 0) {

      ### Es
      Es <- stat$Es_com + t(p_i) %*% X[stat$idxMissingObs, Ys, drop=FALSE]
      ## mapX2Ys <- rep(NA, ncol(X))
      ## mapX2Ys[Ys] <- 1:length(Ys)


      ### Ess
      Ess <- stat$Ess_com + Reduce("+", lapply(1:nrow(levels_Is),
                                               function(k, idxMissingObs, mapAllObs2MissingObs)
                                                 Reduce("+" , lapply(idxMissingObs,
                                                                     function(i, k, mapAllObs2MissingObs) {
                                                                       p_i[mapAllObs2MissingObs[i], k]*t(X[i, Ys, drop=FALSE]) %*% X[i, Ys, drop=FALSE]
                                                                     }, k, mapAllObs2MissingObs)
                                                       ),
                                               stat$idxMissingObs, stat$mapAllObs2MissingObs)
                                  )

      ### ssd
      ssd <- Ess - t(Es)%*%(Es/E_n)

      for (k in 1:nrow(levels_Is)) {
        j <- index[[k]]
        ### WATCH OUT: Es[k, Ys] REPLACED BY Es[k, ] SINCE LOOKS LIKE ALL Ys ARE USED ALONG THE COLUMNS FROM Es
        bar_y[j, Ys_Sigma] <- matrix(rep(Es[k, ]/E_n[k], length(j)), byrow=TRUE, ncol=length(Ys))##, dimnames=list(j, Ys))
      }
      h[ , Ys_Sigma] <- t(n*solve(ssd)%*%t(bar_y[ ,Ys_Sigma]))
      K[Ys_Sigma, Ys_Sigma] <- n*solve(ssd)
    } else {  ### Is != emptyset & Ys == emptyset
      ssd <- 1
    }

  } else {
    ### Is == emptyset & Ys == emptyset
    if (length(Ys)==0) {
      m <- 1
      ssd <- 1
    } else { ### Is == emptyset & Ys != emptyset
      m <- n
      ssd <- (n - 1)*cov(X[ , Ys, drop = FALSE])
      h[, Ys_Sigma] <- matrix(rep(n*solve(ssd) %*% colMeans(X[, Ys, drop=FALSE]), nrow(levels_I)), ncol=length(Ys), byrow=TRUE)
      K[Ys_Sigma, Ys_Sigma] <- n*solve(ssd)
    }
  }

  list(ssd=ssd, K=K, h=h, m=m, bar_y=bar_y)
}

## convergence: calculate convergence diagnostic of the EM algorithm by comparing
## the updated values of the moment parameteres to the moment parameters of the
## previous iteration
convergence <- function(Sigma_update, mu_update, m_update, Sigma, mu, m) {

  delta_m <- abs(m_update - m)/sqrt(m_update + 1)
  delta_mu <- abs(mu_update - mu)/sqrt(diag(Sigma_update))
  delta_Sigma <- abs(Sigma_update - Sigma)
  for (i in 1:ncol(Sigma)) {
    for (j in 1:nrow(Sigma)) {
      delta_Sigma[i,j] <- delta_Sigma[i,j]/sqrt(Sigma_update[i,i]*Sigma_update[j,j] + Sigma_update[i,j]^2)
    }
  }

  list(mdiff=max(max(delta_m), max(delta_mu), max(delta_Sigma)), m=m_update, mu=mu_update, Sigma=Sigma_update)
}

## calculate ssd matrices for H0 and H1 using the EM algorithm
.ssdStatsEM <- function(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
                        I, mapX2I, nLevels, Y, mapX2Y, i, j, Q, tol=0.01) {

  if (length(I) > 0) {
    levels_I <- expand.grid(sapply(nLevels[I], seq_len, simplify=FALSE))
    dimnames(levels_I) <- list(1:nrow(levels_I), I)
  } else {
    levels_I <- matrix(1, ncol=1, nrow=1, dimnames=list(1,1))
  }

  p0 <- rep(1/nrow(levels_I), nrow(levels_I))
  mu0 <- matrix(0, ncol=length(Y), nrow=nrow(levels_I))
  ## colnames(mu0) <- Y
  Sigma0 <- diag(length(Y))
  ## dimnames(Sigma0) <- list(Y, Y)

  p <- p0
  mu <- mu0
  Sigma <- Sigma0
  n <- nrow(X)

  I_j <- intersect(I, c(i, Q))
  Y_j <- intersect(Y, c(i, Q))
  levels_I_j <- comStat_j <- c()
  if (length(I_j) > 0) {
    levels_I_j <- as.matrix(unique(levels_I[, mapX2I[I_j], drop=FALSE]))
    ## colnames(levels_I_j) <- I_j
    comStat_j <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
                          Is=I_j, Ys=Y_j, levels_Is=levels_I_j)
    if (any(comStat_j$n_com == 0))
      stop("Some joint levels of the discrete variables involved lack observations (n(i) == 0).")

    ## cat("Es_com_j=\n")
    ## print(comStat_j$Es_com)
    ## cat("Ess_com_j=\n")
    ## print(comStat_j$Ess_com)
    ## cat("n_com_j=", comStat_j$n_com,"\n")
  }

  I_i <- intersect(I, c(j, Q))
  Y_i <- intersect(Y, c(j, Q))
  levels_I_i <- comStat_i <- c()

  if (length(I_i) > 0) {
    levels_I_i <- as.matrix(unique(levels_I[, mapX2I[I_i], drop=FALSE]))
    ## colnames(levels_I_i) <- I_i
    comStat_i <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
                          Is=I_i, Ys=Y_i, levels_Is=levels_I_i)
    if (any(comStat_i$n_com == 0))
      stop("Some joint levels of the discrete variables involved lack observations (n(i) == 0).")

    ## cat("Es_com_i=\n")
    ## print(comStat_i$Es_com)
    ## cat("Ess_com_i=\n")
    ## print(comStat_i$Ess_com)
    ## cat("n_com_i=", comStat_i$n_com,"\n")
  }

  I_ij <- intersect(I, Q)
  Y_ij <- intersect(Y, Q)
  levels_I_ij <- comStat_ij <- c()
  if (length(I_ij) > 0) {
    levels_I_ij <- as.matrix(unique(levels_I[, mapX2I[I_ij], drop=FALSE]))
    ## colnames(levels_I_ij) <- I_ij
    comStat_ij <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
                           Is=I_ij, Ys=Y_ij, levels_Is=levels_I_ij)
    if (any(comStat_ij$n_com == 0))
      stop("Some joint levels of the discrete variables involved lack observations (n(i) == 0).")

    ## cat("Es_com_ij=\n")
    ## print(comStat_ij$Es_com)
    ## cat("Ess_com_ij=\n")
    ## print(comStat_ij$Ess_com)
    ## cat("n_com_ij=", comStat_ij$n_com,"\n")
  }

  mdiff <- 1
  while (mdiff > tol) {
    sufstat_j <- stat_mis(X, I_j, Y_j, levels_I_j, comStat_j, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma)
    sufstat_i <- stat_mis(X, I_i, Y_i, levels_I_i, comStat_i, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma)
    sufstat_ij <- stat_mis(X, I_ij, Y_ij, levels_I_ij, comStat_ij, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma)

    K <- sufstat_j$K + sufstat_i$K - sufstat_ij$K
    h <- sufstat_j$h + sufstat_i$h - sufstat_ij$h
    m <- sufstat_j$m*sufstat_i$m/sufstat_ij$m

    if (all(K == matrix(0,ncol=length(Y), nrow=length(Y)))){
      Sigma_update <- 0
    } else {
      Sigma_update <- solve(K)
    }
    mu_update <- t(Sigma_update%*%t(h))
    conv <- convergence(Sigma_update=Sigma_update, mu_update=mu_update, m_update=m, Sigma=Sigma, mu=mu, m=p*n)
    mdiff <- conv$mdiff
    p <- conv$m/n
    mu <- conv$mu
    Sigma <- conv$Sigma
  }

  mdiff <- 1
  p <- p0
  mu <- mu0
  Sigma <- Sigma0
  comStat <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs, Is=I, Ys=Y, levels_Is=levels_I)

  while (mdiff > tol) {
    sufstat <- stat_mis(X, I, Y, levels_I, comStat, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma)

    ssd <- sufstat$ssd
    bar_y <- sufstat$bar_y
    m <- sufstat$m

    conv <- convergence(Sigma_update=ssd/n, mu_update=bar_y, m_update=m, Sigma=Sigma, mu=mu, m=p*n)
    mdiff <- conv$mdiff
    p <- conv$m/n
    mu <- conv$mu
    Sigma <- conv$Sigma
  }

  list(ssd_j=sufstat_j$ssd, ssd_i=sufstat_i$ssd, ssd_ij=sufstat_ij$ssd, ssd=sufstat$ssd)
}

.qpCItestHMGM <- function(X, I, nLevels, Y, ssdMat, mapX2ssdMat, i, j, Q,
                          exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                          R.code.only=FALSE ) {

  if (!is.null(ssdMat)) {
    p <- (d <- dim(ssdMat))[1]
    if (p != d[2] || !isSymmetric(ssdMat))
      stop("ssdMat is not squared and symmetric. Is it really an ssd matrix?\n")
    if (class(ssdMat) != "SsdMatrix")
      stop("qpgraph:::.qpCItestHMGM: the ssdMat argument should be an object of class SsdMatrix\n")
  }

  if (all(!is.na(match(c(i,j), I))))
    stop("i and j cannot be both discrete at the moment")

  if (!R.code.only) {
    return(qpgraph:::.qpFastCItestHMGM(X, I, nLevels, Y, ssdMat, mapX2ssdMat,
                                       i, j, Q, exact.test, use, tol))
  }

  if (!is.na(match(j, I))) { ## if any of (i,j) is discrete, it should be i
    tmp <- i
    i <- j
    j <- tmp
  }
  ## cat ("\n", i, "ci", j,"\n")

  I <- intersect(I, c(i, Q))
  Y <- intersect(Y, c(i, j, Q))

  ssd <- ssd_i <- ssd_j <- ssd_ij <- diag(2) 
  n_co <- n <- nrow(X)

  ## build logical mask of missing observations
  missingMask <- apply(X[, c(I, Y), drop=FALSE], 1, function(x) any(is.na(x)))
  missingData <- any(missingMask)

  if (!missingData || use == "complete.obs") { ## either no missing data or use complete obs
    n_co <- n - sum(missingMask)

    if (!missingData && !is.null(ssdMat) && length(I) == 0)
      ssd <- ssdMat[mapX2ssdMat[Y], mapX2ssdMat[Y], drop=FALSE]
    else
      ssd <- qpgraph:::.ssdStatsCompleteObs(X, I, Y, missingMask)
    ## cat("ssd:\n")
    ## print(as.matrix(ssd))

    ## ssd_i = ssd_Gamma when i is discrete or ssd_{Gamma\i} when i is continuous
    if (!missingData && !is.null(ssdMat) && length(setdiff(I, i)) == 0)
      ssd_i <- ssdMat[mapX2ssdMat[setdiff(Y, i)], mapX2ssdMat[setdiff(Y, i)], drop=FALSE]
    else
      ssd_i <- qpgraph:::.ssdStatsCompleteObs(X, setdiff(I, i), setdiff(Y, i), missingMask)
    ## cat("ssd_i:\n")
    ## print(as.matrix(ssd_i))

    if (length(setdiff(Y, j)) > 0) {
      ssd_j <- qpgraph:::.ssdStatsCompleteObs(X, I, setdiff(Y, j), missingMask)
      ## cat("ssd_j:\n")
      ## print(ssd_j)
      if (length(setdiff(Y, c(i,j))) > 0) {
        if (!missingData && !is.null(ssdMat) && length(setdiff(I, i)) == 0)
          ssd_ij <- ssdMat[mapX2ssdMat[setdiff(Y, c(i, j))],
                           mapX2ssdMat[setdiff(Y, c(i, j))], drop=FALSE]
        else
          ssd_ij <- qpgraph:::.ssdStatsCompleteObs(X, setdiff(I, i), setdiff(Y, c(i, j)), missingMask)
        ## cat("ssd_j:\n")
        ## print(ssd_ij)
      }
    }
  } else { ## missing data and should use the EM algorithm
    missingMask2 <- apply(X[, Y, drop=FALSE], 1, function(x) any(is.na(x)))
    if (length(I) == 0 || any(missingMask2))
      stop("EM not implemented yet for missing values in continuous variables. Please set use=\"complete.obs\"\n")

    mapX2Y <- rep(NA, ncol(X))
    mapX2Y[Y] <- 1:length(Y)
    mapX2I <- rep(NA, ncol(X))
    mapX2I[I] <- 1:length(I)

    idxMissingObs <- which(missingMask)
    idxCompleteObs <- setdiff(1:n, idxMissingObs)
    mapAllObs2MissingObs <- rep(NA, n)
    mapAllObs2MissingObs[idxMissingObs] <- 1:length(idxMissingObs)
    ssdMats <- qpgraph:::.ssdStatsEM(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
                                       I, mapX2I, nLevels, Y, mapX2Y, i, j, Q, tol)
    ssd <- as.matrix(ssdMats$ssd)
    ssd_i <- as.matrix(ssdMats$ssd_i)
    ssd_j <- as.matrix(ssdMats$ssd_j)
    ssd_ij <- as.matrix(ssdMats$ssd_ij)
  }
 
  ssd <- determinant(ssd)        ## WATCH OUT! when using Matrix::determinant(..., logarithm=TRUE)
  ssd_i <- determinant(ssd_i)    ## $modulus is always 0, don't know why. since this is its default
  ssd_j <- determinant(ssd_j)    ## this argument is not being put explicitly in the call
  ssd_ij <- determinant(ssd_ij)  ## keep an eye in case the default ever changes
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
      param <- c(a=a, b=b, n=ifelse(use == "em", n, n_co))
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
      param <- c(df=df, n=ifelse(use == "em", n, n_co))
      p.value <- c(greater = 1 - pchisq(lr, df=df, lower.tail=TRUE))
      n.value <- c("-n log Lambda" = 0)
      method  <- "Conditional independence test for homogeneous mixed data using an asymptotic likelihood ratio test"
      alt <- "greater"
    }
  }

  RVAL <- list(statistic=stat,
               parameter=param,
               p.value=if (use != "em") p.value else NA_real_, ## p-values currently not valid with EM
               estimate=NULL,
               null.value=n.value,
               alternative=alt,
               method=method,
               data.name=sprintf("%s and %s given {%s}", names(i), names(j), paste(names(Q), collapse=", ")))
  class(RVAL) <- "htest"

  RVAL
}

setGeneric("qpAllCItests", function(X, ...) standardGeneric("qpAllCItests"))

setMethod("qpAllCItests", signature(X="matrix"),
          function(X, I=NULL, Q=NULL, pairup.i=NULL, pairup.j=NULL,
                   long.dim.are.variables=TRUE, exact.test=TRUE,
                   use=c("complete.obs", "em"), tol=0.01, return.type=c("p.value", "statn", "all"),
                   verbose=TRUE, R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                   nAdj2estimateTime=10) {

            use <- match.arg(use)
            return.type <-  match.arg(return.type)

            startTime <- c(user.self=0, sys.self=0, elapsed=0, user.child=0, sys.child=0)
            class(startTime) <- "proc_time"
            if (estimateTime)
              startTime <- proc.time()

            if (clusterSize > 1 && R.code.only)
              stop("Using a cluster (clusterSize > 1) only works with R.code.only=FALSE\n")

            if (clusterSize > 1 &&
               (!qpgraph:::.qpIsPackageInstalled("snow") || !qpgraph:::.qpIsPackageInstalled("Rmpi")))
              stop("Using a cluster (clusterSize > 1) requires first installing packages 'snow' and 'Rmpi'\n")

            if (long.dim.are.variables &&
                sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
              X <- t(X)

            if (is.null(colnames(X))) 
              colnames(X) <- 1:ncol(X)

            qpgraph:::.qpAllCItests(X, I, Q, pairup.i, pairup.j,
                                    exact.test, use, tol, return.type,
                                    verbose, R.code.only, clusterSize,
                                    startTime, nAdj2estimateTime)
          })

.qpAllCItests <- function(X, I=NULL, Q=NULL, pairup.i=NULL, pairup.j=NULL,
                          exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                          return.type=c("p.value", "statn", "all"), verbose=TRUE,
                          R.code.only=FALSE, clusterSize=1, startTime, nAdj2estimateTime=10) {

  cl <- NULL

  if (use == "em" && !R.code.only)
    stop("use=\"em\" does not work yet with R.code.only=FALSE\n")
 
  if (class(clusterSize)[1] == "numeric" || class(clusterSize)[1] == "integer") {
    if (clusterSize > 1) {
      ## copying ShortRead's strategy, 'get()' are to quieten R CMD check, and for no other reason
      ## makeCl <- get("makeCluster", mode="function")
      ## clSetupRNG <- get("clusterSetupRNG", mode="function")
      ## clEvalQ <- get("clusterEvalQ", mode="function")
      ## clExport <- get("clusterExport", mode="function")
      ## clApply <- get("clusterApply", mode="function")
      ## stopCl <- get("stopCluster", mode="function")
      ## clCall <- get("clusterCall", mode="function")
      ## clOpt <- get("getClusterOption", mode="function")

      if (startTime["elapsed"] == 0)
        message("Testing conditional independences using a ", snow::getClusterOption("type"),
                " cluster of ", clusterSize, " nodes\n")
      else
        message("Estimating time of testing conditional independences using a ", snow::getClusterOption("type"),
                " cluster of ", clusterSize, " nodes\n")

      ## REPLACEMENT OF PACKAGE rlecuyer
      ## cl <- makeCl(clusterSize, snowlib=system.file(package="qpgraph"))
      ## clSetupRNG(cl)
      cl <- parallel::makeCluster(clusterSize, type="MPI", snowlib=system.file(package="qpgraph"))
      parallel::clusterSetRNGStream(cl)
      res <- parallel::clusterEvalQ(cl, require(qpgraph, quietly=TRUE))
      if (!all(unlist(res))) {
        parallel::stopCluster(cl)
        stop("The package 'qpgraph' could not be loaded in some of the nodes of the cluster")
      }
      assign("clusterSize", clusterSize, envir=.GlobalEnv)
      parallel::clusterExport(cl, list("clusterSize"))
      rm("clusterSize", envir=.GlobalEnv)
      parallel::clusterApply(cl, 1:clusterSize, function(x) assign("clusterRank", x, envir=.GlobalEnv))
    }
  } else {
    if (!is.na(match("cluster", class(clusterSize))))
      cl <- clusterSize
    else
      stop("Unknown class for argument clusterSize:", class(clusterSize))
  }

  ## X the matrix of data with columns as r.v. and rows as multivariate observations

  var.names <- colnames(X)
  n.var <- ncol(X)
  N <- nrow(X)

  ## check that the q parameter has proper values

  q <- length(Q)

  if (q > n.var - 2)
    stop(paste("q=",q," > p-2=", n.var-2))

  if (q < 0)
    stop(paste("q=",q," < 0"))

  if (q > N - 3)
    stop(paste("q=",q," > n-3=",N-3))

  ## check whether there are discrete variables and whether they're properly set

  nLevels <- rep(NA_integer_, times=ncol(X))
  Y <- NULL
  if (!is.null(I)) {
    if (is.character(I)) {
      if (any(is.na(match(I, var.names))))
        stop("Some variables in I do not form part of the variable names of the data in X\n")
      I <- match(I, var.names)
    } else {
      if (any(is.na(match(I, 1:n.var))))
        stop("Some variables in I do not form part of the variables of the data in X\n")
    }
    nLevels[I] <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
    if (any(nLevels[I] == 1))
      stop(sprintf("Discrete variable %s has only one level", colnames(X)[I[nLevels[I]==1]]))

    Y <- (1:n.var)[-I]
  }

  if ((!is.null(pairup.i) && is.null(pairup.j)) ||
      (is.null(pairup.i) && !is.null(pairup.j)))
    stop("pairup.i and pairup.j should both either be set to NULL or contain subsets of variables\n")

  if (is.null(pairup.i))
    pairup.i <- 1:n.var
  else {
    if (is.character(pairup.i)) {
      if (any(is.na(match(pairup.i, var.names))))
        stop("Some variables in pairup.i do not form part of the variable names of the data in X\n")
      pairup.i <- match(pairup.i, var.names)
    }
  }

  if (is.null(pairup.j)) {
    pairup.j <- 1:n.var
    if (!is.null(I)) { ## by now, interactions between discrete variables are not considered
        pairup.j <- (1:n.var)[-I]
    }
  } else {
    if (is.character(pairup.j)) {
      if (any(is.na(match(pairup.j, var.names))))
        stop("Some variables in pairup.j do not form part of the variable names of the data in X\n")
      pairup.j <- match(pairup.j, var.names)
    }
  }

  if (!is.null(Q)) {
    if (is.character(Q)) {
      if (any(is.na(match(Q, var.names))))
        stop("Some variables in Q do not form part of the variable names of the data\n")
      Q <- match(Q, var.names)
    } else {
      if (any(is.na(match(Q, 1:n.var))))
        stop("Some variables in Q do not form part of the variables of the data\n")
    }

    ## variables in Q are removed from the pairs for which CI tests are performed
    pairup.i <- setdiff(pairup.i, Q)
    pairup.j <- setdiff(pairup.j, Q)
  }

  ## pair the two sets pairup.i and pairup.j without pairing the same variable
  l.pairup.i <- length(pairup.i)
  l.pairup.j <- length(pairup.j)
  l.int <- length(intersect(pairup.i, pairup.j))
  l.pairup.i.noint <- l.pairup.i - l.int
  l.pairup.j.noint <- l.pairup.j - l.int
  n.adj <- l.int * l.pairup.j.noint + l.int * l.pairup.i.noint +
           l.pairup.i.noint * l.pairup.j.noint + l.int * (l.int - 1) / 2

  pairup.ij.int <- intersect(pairup.i, pairup.j)
  pairup.i.noint <- setdiff(pairup.i, pairup.ij.int)
  pairup.j.noint <- setdiff(pairup.j, pairup.ij.int)

  pvstatn <- list(p.value=NA, statistic=NA, n=NA)

  if (!R.code.only) {
    elapsedTime <- 0
    if (startTime["elapsed"] > 0) {
      elapsedTime <- (proc.time() - startTime)["elapsed"]
      startTime <- proc.time()
    }

    if (is.null(cl)) { ## single-processor execution

      cit <- qpgraph:::.qpFastAllCItests(X, I, Y, Q, pairup.i.noint,
                                         pairup.j.noint, pairup.ij.int,
                                         exact.test, use, tol, return.type, verbose,
                                         startTime["elapsed"], nAdj2estimateTime)

      if (startTime["elapsed"] == 0) {
        if (return.type == "all" || return.type == "p.value") {
          pvstatn$p.value <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                                 Dimnames=list(var.names, var.names), x=cit$p.value)
        }
        if (return.type == "all" || return.type == "statn") {
          pvstatn$statistic <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                                   Dimnames=list(var.names, var.names), x=cit$statistic)
          pvstatn$n <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                           Dimnames=list(var.names, var.names), x=as.numeric(cit$n))
        }
      } else ## completion time estimation comes directly as a vector
        pvstatn <- cit

    } else {           ## use a cluster !
      ## clCall <- get("clusterCall", mode="function")
      valIdx <- list()
      if (verbose && startTime["elapsed"] == 0) { ## no cluster progress-call when only estimating time
        valIdx <- clPrCall(cl, qpgraph:::.qpFastAllCItestsPar, n.adj, X, I, Y, Q,
                           pairup.i.noint, pairup.j.noint, pairup.ij.int,
                           exact.test, use, tol, return.type, verbose, FALSE, nAdj2estimateTime)
      } else {
        valIdx <- parallel::clusterCall(cl, qpgraph:::.qpFastAllCItestsPar, X, I, Y, Q,
                                        pairup.i.noint, pairup.j.noint, pairup.ij.int,
                                        exact.test, use, tol, return.type, verbose, startTime["elapsed"] > 0,
                                        nAdj2estimateTime)
      }

      if (startTime["elapsed"] > 0) {
        ## the following calculation makes important part of the estimation of the time
        ## it assumes that the estimated time per processor is stored on the first position of 'p.value'
        ## and uses the median of the times estimated for each processor to try to be robust against
        ## fluctuations on the execution time taken in some processors
        elapsedTime <- elapsedTime + median(sapply(valIdx, function(x) x$p.value[1]))
        startTime <- proc.time()
      }

      if (class(clusterSize)[1] == "numeric" || class(clusterSize)[1] == "integer")
        parallel::stopCluster(cl)

        if (return.type == "all" || return.type == "p.value") {
          pvstatn$p.value <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                                 Dimnames=list(var.names, var.names),
                                 x=rep(as.double(NA), n.var*(n.var-1)/2+n.var)) 
          pvstatn$p.value@x[do.call("c", lapply(valIdx, function(x) x$idx))] <-
            do.call("c", lapply(valIdx, function(x) x$p.value))
        }

        if (return.type == "all" || return.type == "statn") {
          pvstatn$statistic <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                                   Dimnames=list(var.names, var.names),
                                   x=rep(as.double(NA), n.var*(n.var-1)/2+n.var)) 
          pvstatn$statistic@x[do.call("c", lapply(valIdx, function(x) x$idx))] <-
            do.call("c", lapply(valIdx, function(x) x$statistic))

          pvstatn$n <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                           Dimnames=list(var.names, var.names),
                           x=rep(as.double(NA), n.var*(n.var-1)/2+n.var)) 
          pvstatn$n@x[do.call("c", lapply(valIdx, function(x) x$idx))] <-
            do.call("c", lapply(valIdx, function(x) as.numeric(x$n)))
        }

      if (startTime["elapsed"] > 0) {
        elapsedTime <- elapsedTime + (proc.time() - startTime)["elapsed"]
        d <- as.vector(floor(elapsedTime / (24*3600)))
        h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
        m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
        s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
        pvstatn <- c(days=d, hours=h, minutes=m, seconds=s)
      }
    }

    return(pvstatn)
  }

  missingMask <- apply(X[, , drop=FALSE], 1, function(x) any(is.na(x)))
  missingData <- any(is.na(missingMask))

  S <- ssd <- mapX2ssd <- NULL
  if (!missingData) {
    if (!is.null(I)) {  ## calculate the uncorrected sum of squares and deviations
      ssd <- qpCov(X[, Y, drop=FALSE], corrected=FALSE)
      mapX2ssd <- match(var.names, colnames(ssd))
      ## names(mapX2ssd) <- colnames(X) ## is this necessary?
    } else             ## calculate sample covariance matrix
      S <- qpCov(X)
  }

  ## return an efficiently stored symmetric matrix
  if (return.type == "all" || return.type == "p.value") {
    pvstatn$p.value <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                           Dimnames=list(var.names, var.names),
                           x=rep(as.double(NA), n.var*(n.var-1)/2+n.var))
  }

  if (return.type == "all" || return.type == "statn") {
    pvstatn$statistic <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                             Dimnames=list(var.names, var.names),
                             x=rep(as.double(NA), n.var*(n.var-1)/2+n.var))
    pvstatn$n <- new("dspMatrix", Dim=as.integer(c(n.var, n.var)),
                     Dimnames=list(var.names, var.names),
                     x=rep(as.double(NA), n.var*(n.var-1)/2+n.var))
  }

  elapsedTime <- 0
  if (startTime["elapsed"] > 0) {
    elapsedTime <- (proc.time() - startTime)["elapsed"]
    startTime <- proc.time()
  }

  ppct <- -1
  k <- 0
  pb <- NULL
  if (verbose && elapsedTime == 0)
    pb <- txtProgressBar(style=3)

  cit <- NA ### build return object depending on return.type CONTINUE HERE !!!!

  ## intersection variables against ij-exclusive variables
  for (i in pairup.ij.int) {
    for (j in c(pairup.i.noint,pairup.j.noint)) {

      if (is.null(I)) {
        Xsub <- X[, c(i, j, Q), drop=FALSE]
        S <- qpCov(Xsub)
        cit <- qpgraph:::.qpCItest(S, 1L, 2L, 2L+seq(along=Q), R.code.only=TRUE)
      } else
        cit <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                       exact.test, use, tol, R.code.only=TRUE)

      if (return.type == "all" || return.type == "p.value")
        pvstatn$p.value[j,i] <- pvstatn$p.value[i,j] <- cit$p.value
      if (return.type == "all" || return.type == "statn") {
        pvstatn$statistic[j,i] <- pvstatn$statistic[i,j] <- cit$statistic
        pvstatn$n[j,i] <- pvstatn$n[i,j] <- cit$param["n"]
      }

      k <- k + 1
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
      pct <- floor((k * 100) / n.adj)
      if (pct != ppct && verbose && elapsedTime == 0) {
        setTxtProgressBar(pb, pct/100)
        ppct <- pct
      }
    }
    if (elapsedTime > 0 && k == nAdj2estimateTime)
      break;
  }

  ## i-exclusive variables against j-exclusive variables
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in pairup.i.noint) {
      for (j in pairup.j.noint) {

        if (is.null(I)) {
          Xsub <- X[, c(i, j, Q), drop=FALSE]
          S <- qpCov(Xsub)
          cit <- qpgraph:::.qpCItest(S, 1L, 2L, 2L+seq(along=Q), R.code.only=TRUE)
        } else
          cit <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                                         exact.test, use, tol, R.code.only=TRUE)

        if (return.type == "all" || return.type == "p.value")
          pvstatn$p.value[j,i] <- pvstatn$p.value[i,j] <- cit$p.value
        if (return.type == "all" || return.type == "statn") {
          pvstatn$statistic[j,i] <- pvstatn$statistic[i,j] <- cit$statistic
          pvstatn$n[j,i] <- pvstatn$n[i,j] <- cit$param["n"]
        }

        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  ## intersection variables against themselves (avoiding pairing of the same)
  if (elapsedTime == 0 || k < nAdj2estimateTime) {
    for (i in seq(along=pairup.ij.int[-1])) {
      i2 <- pairup.ij.int[i]

      for (j in (i+1):l.int) {
        j2 <- pairup.ij.int[j]

        if (is.null(I)) {
          Xsub <- X[, c(i, j, Q), drop=FALSE]
          S <- qpCov(Xsub)
          cit <- qpgraph:::.qpCItest(S, 1L, 2L, 2L+seq(along=Q), R.code.only=TRUE)
        } else
          cit <- qpgraph:::.qpCItestHMGM(X, I, nLevels, Y, ssd, mapX2ssd, i2, j2, Q,
                                         exact.test, use, tol, R.code.only=TRUE)

        if (return.type == "all" || return.type == "p.value")
          pvstatn$p.value[j2,i2] <- pvstatn$p.value[i2,j2] <- cit$p.value
        if (return.type == "all" || return.type == "statn") {
          pvstatn$statistic[j2,i2] <- pvstatn$statistic[i2,j2] <- cit$statistic
          pvstatn$n[j2,i2] <- pvstatn$n[i2,j2] <- cit$param["n"]
        }

        k <- k + 1
        if (elapsedTime > 0 && k == nAdj2estimateTime)
          break;
        pct <- floor((k * 100) / n.adj)
        if (pct != ppct && verbose && elapsedTime == 0) {
          setTxtProgressBar(pb, pct/100)
          ppct <- pct
        }
      }
      if (elapsedTime > 0 && k == nAdj2estimateTime)
        break;
    }
  }

  if (verbose && elapsedTime == 0) {
    close(pb)
  }

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + ((proc.time()-startTime)["elapsed"]/k) * n.adj
    startTime <- proc.time()
  }

  ## this is necessary till we find out how to properly assign values in a dspMatrix
  ## nrrMatrix <- as(nrrMatrix, "dspMatrix")

  if (elapsedTime > 0) {
    elapsedTime <- elapsedTime + (proc.time()-startTime)["elapsed"]
    d <- as.vector(floor(elapsedTime / (24*3600)))
    h <- as.vector(floor((elapsedTime-d*24*3600)/3600))
    m <- as.vector(floor((elapsedTime-d*24*3600-h*3600)/60))
    s <- as.vector(ceiling(elapsedTime-d*24*3600-h*3600-m*60))
    pvstatn <- c(days=d, hours=h, minutes=m, seconds=s)
  }

  return(pvstatn)
}


.qpFastCItestStd <- function(S, i, j, Q) {
  return(.Call("qp_fast_ci_test_std", S@ssd@x, nrow(S), as.integer(S@n), i, j, Q))
}

.qpFastCItestOpt <- function(S, n, i, j, Q) {
  return(.Call("qp_fast_ci_test_opt", S@ssd@x, nrow(S), as.integer(S@n), i, j, Q))
}

.qpFastCItestHMGM <- function(X, I, nLevels, Y, ssd, mapX2ssd, i, j, Q,
                              exact.test, use, tol) {
  x <- NULL
  if (!is.null(ssd)) {
    if (class(ssd) != "SsdMatrix")
      stop("qpgraph:::.qpFastCItestHMGM: the ssd argument should be an object of class SsdMatrix\n")
    x <- ssd@ssd@x
  }

  return(.Call("qp_fast_ci_test_hmgm", X, I, nLevels, Y, x, as.integer(mapX2ssd),
               i, j, Q, as.integer(exact.test),
               as.integer(factor(use, levels=c("complete.obs", "em"))), tol))
}

.qpFastAllCItests <- function(X, I, Y, Q, pairup.i.noint, pairup.j.noint,
                              pairup.ij.int, exact.test, use, tol, return.type,
                              verbose, startTime, nAdj2estimateTime) {
  nLevels <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))
  return(.Call("qp_fast_all_ci_tests", X, as.integer(I), as.integer(nLevels),
                                      as.integer(Y), as.integer(Q),
                                      as.integer(pairup.i.noint), as.integer(pairup.j.noint),
                                      as.integer(pairup.ij.int), as.integer(exact.test),
                                      as.integer(factor(use, levels=c("complete.obs", "em"))), tol,
                                      as.integer(factor(return.type, levels=c("p.value", "statn", "all"))),
                                      as.integer(verbose), as.double(startTime),
                                      as.integer(nAdj2estimateTime), .GlobalEnv))
}

.qpFastAllCItestsPar <- function(X, I, Y, Q, pairup.i.noint, pairup.j.noint,
                                 pairup.ij.int, exact.test, use, tol, return.type,
                                 verbose, estimateTime, nAdj2estimateTime) {
  ## clOpt <- get("getClusterOption", mode="function")
  myMaster <- snow::getClusterOption("masterNode")

  startTime <- 0
  if (estimateTime)
    startTime <- proc.time()["elapsed"]

  nLevels <- apply(X[, I, drop=FALSE], 2, function(x) nlevels(as.factor(x)))

  ## clusterRank and clusterSize should have been defined by the master node
  return(.Call("qp_fast_all_ci_tests_par", X, as.integer(I), as.integer(nLevels),
                                           as.integer(Y), as.integer(Q),
                                           as.integer(pairup.i.noint), as.integer(pairup.j.noint),
                                           as.integer(pairup.ij.int), as.integer(exact.test),
                                           as.integer(factor(use, levels=c("complete.obs", "em"))), tol,
                                           as.integer(factor(return.type, levels=c("p.value", "statn", "all"))),
                                           as.integer(verbose), as.double(startTime),
                                           as.integer(nAdj2estimateTime), as.integer(get("clusterRank")),
                                           as.integer(get("clusterSize")), myMaster, .GlobalEnv))
}
