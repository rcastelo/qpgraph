## virtual class graphParam to hold parameters to simulate graphs
setClass("graphParam",
         representation(p="integer",
                        labels="character"),
         prototype(p=5L,
                   labels=as.character(1:5)))

## virtual class markedGraphParam to hold parameters to simulate marked graphs
## where a marked graph is a graph with marked vertices
setClass("markedGraphParam",
         representation(pI="integer",
                        pY="integer",
                        Ilabels="character",
                        Ylabels="character"),
         prototype(pI=1L,
                   pY=4L,
                   Ilabels="I1",
                   Ylabels=paste0("Y", 2:5)))

## class erGraphParam to hold parameters to simulate Erdos-Renyi graphs
setClass("erGraphParam",
         contains="graphParam",
         representation(m="integer",
                        prob="numeric"),
         prototype(p=5L,
                   m=5L,
                   prob=NA_real_,
                   labels=as.character(1:5)))

## class dRegularGraphParam to hold parameters to simulate d-regular graphs
setClass("dRegularGraphParam",
         contains="graphParam",
         representation(d="integer",
                        exclude="integer"),
         prototype(p=5L,
                   d=2L,
                   exclude=as.integer(NULL),
                   labels=as.character(1:5)))

## class erMarkedGraphParam to hold parameters to simulate Erdos-Renyi marked graphs
setClass("erMarkedGraphParam",
         contains=c("markedGraphParam", "erGraphParam"),
         representation(),
         prototype(p=5L,
                   pI=1L,
                   pY=4L,
                   Ilabels="I1",
                   Ylabels=paste0("Y", 2:5),
                   labels=c("I1", paste0("Y", 2:5)),
                   m=5L,
                   prob=NA_real_))

## class dRegularMarkedGraphParam to hold parameters to simulate d-regular marked graphs
setClass("dRegularMarkedGraphParam",
         contains=c("markedGraphParam", "dRegularGraphParam"),
         representation(),
         prototype(p=5L,
                   pI=1L,
                   pY=4L,
                   Ilabels="I1",
                   Ylabels=paste0("Y", 2:5),
                   labels=c("I1", paste0("Y", 2:5)),
                   d=2L,
                   exclude=as.integer(NULL)))

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
                   sigma=as(as(as(diag(1:5), "dMatrix"), "symmetricMatrix"), "packedMatrix")))

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
                   sigma=as(as(as(diag(1:5), "dMatrix"), "symmetricMatrix"), "packedMatrix"),
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


setOldClass(c("bc", "cross"))
setOldClass("map")

## class eQTLcrossParam to hold parameters to simulate eQTLcross objects
setClass("eQTLcrossParam",
         representation(map="map",
                        type="character",
                        cis="numeric",
                        trans="integer",
                        cisr="numeric",
                        d2m="numeric",
                        networkParam="graphParam"))

## class eQTLcross to hold an experimental cross involving genotype markers
## and gene expression profiles with some underlying expression
## quantitative trait loci (eQTL) and some underlying regulatory network between genes
setClass("eQTLcross",
         representation(map="map",
                        genes="matrix",
                        model="HMgmm",
                        type="character"))

setClass("eQTLnetworkEstimationParam",
         representation(ggData="matrix",
                        geneticMap="ANY",
                        physicalMap="ANY",
                        organism="character",
                        genome="Seqinfo",
                        geneAnnotation="GRanges",
                        geneAnnotationTable="character",
                        dVars="character"))

setClass("qpGraph",
         representation(p="integer",
                        q="integer",
                        n="integer",
                        epsilon="numeric",
                        g="graphBAM"))

setClass("eQTLnetwork",
         representation(geneticMap="ANY",
                        physicalMap="ANY",
                        organism="character",
                        genome="Seqinfo",
                        geneAnnotation="GRanges",
                        geneAnnotationTable="character",
                        dVars="character",
                        pvaluesG0="dspMatrix",
                        nrr="dspMatrix",
                        modelFormula="formula",
                        rhs="list",
                        qOrders="integer",
                        p.value="numeric",
                        adjustMethod="character",
                        epsilon="numeric",
                        alpha="numeric",
                        qpg="qpGraph"))


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
