## remove co-linearities in the longer dimension

filterCollinearities <- function(X, soft.filter=FALSE, long.dim.are.variables=TRUE) {
  r <- qpPCC(X, long.dim.are.variables=long.dim.are.variables)$R
  r[upper.tri(r)] <- 0
  diag(r) <- 0
  mask <- apply(r, 2, function(x) any(x > 0.99))

  if (!soft.filter) {
   if (long.dim.are.variables &&
       sort(dim(X),decreasing=TRUE,index.return=TRUE)$ix[1] == 1)
     X <- X[!mask, ]
   else
     X <- X[, !mask]
  } else
    return(mask)

  X
}
