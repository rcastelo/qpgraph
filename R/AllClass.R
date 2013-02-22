## class UGgmm to simulate an undirected Gaussian graphical Markov model
setClass("UGgmm",
         representation(p="integer",
                        g="graphBAM",
                        mean="numeric",
                        sigma="dspMatrix"),
         prototype(p=5L,
                   g=graphBAM(as.data.frame(matrix(NA, nrow=0, ncol=3,
                                                   dimnames=list(NULL, c("from", "to", "weight")))),
                              nodes=sprintf("X%d", 1:5)),
                   mean=do.call("names<-", list(rep(0, 5), sprintf("X%d", 1:5))),
                   sigma=as(diag(1:5), "dspMatrix")))

setClass("UGgmmSummary",
         representation(model="UGgmm",
                        density="numeric",
                        degree="integer",
                        macor="numeric",
                        pacor="numeric"))

## class HMgmm to simulate an undirected Gaussian graphical Markov model
setClass("HMgmm",
         representation(pI="integer",
                        pY="integer",
                        g="graphBAM",
                        vtype="factor",
                        dLevels="integer",
                        a="numeric",
                        rho="numeric",
                        sigma="dspMatrix",
                        mean="environment",
                        eta2="environment"),
         prototype(pI=1L,
                   pY=4L,
                   vtype=factor(c("discrete", rep("continuous", 4))),
                   g={g<- graphBAM(as.data.frame(matrix(NA, nrow=0, ncol=3,
                                                  dimnames=list(NULL, c("from", "to", "weight")))),
                                   nodes=c("I01", sprintf("%d", 1:4)))
                      nodeDataDefaults(g, "type") <- "continuous"
                      nodeData(g, "I01", "type") <- "discrete"
                      g},
                   dLevels=2L,
                   a=do.call("names<-", list(rep(0, 4), sprintf("Y%d", 1:4))),
                   rho=0.5,
                   sigma=as(diag(1:5), "dspMatrix"),
                   mean=new.env(parent=emptyenv()),
                   eta2=new.env(parent=emptyenv())))
                   ## eta2=do.call("names<-", list(rep(NA, 4), sprintf("Y%d", 1:4)))))

setClass("HMgmmSummary",
         representation(model="HMgmm",
                        density="numeric",
                        densityIxY="numeric",
                        densityY="numeric",
                        degree="integer",
                        macor="numeric",
                        pacor="numeric",
                        a="numeric"))


## class eQTLcross to hold an experimental cross involving genotype markers
## and gene expression profiles with some underlying expression
## quantitative trait loci (eQTL) and some underlying regulatory network between genes

setOldClass(c("bc", "cross"))
setOldClass("map")
setClass("eQTLcross",
         representation(map="map",
                        genes="matrix",
                        model="HMgmm",
                        type="character"))


## class SsdMatrix to store matrices with sum of squares of deviations (ssd)
## these are dspMatrix objects with an additional slot 'n' indicating the
## sample size from where these ssd matrices were estimated. the main use
## of this extra slot is to inform the user of the sample size when the
## function qpCov() was called with use="complete.obs" on data with missing
## values
setClass("SsdMatrix",
         representation(ssd = "dspMatrix",
                        n = "numeric"),
         prototype(ssd=new("dspMatrix"), n=0),
         contains = "dspMatrix")
