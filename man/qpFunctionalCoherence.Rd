\name{qpFunctionalCoherence}
\alias{qpFunctionalCoherence}
\alias{qpFunctionalCoherence,matrix-method}
\alias{qpFunctionalCoherence,lsCMatrix-method}
\alias{qpFunctionalCoherence,lspMatrix-method}
\alias{qpFunctionalCoherence,lsyMatrix-method}
\alias{qpFunctionalCoherence,list-method}

\title{
Functional coherence estimation
}
\description{
Estimates functional coherence for a given transcriptional regulatory network
specified either as an adjacency matrix with a list of transcription factor gene
identifiers or as a list of transcriptional regulatory modules.
}
\usage{
\S4method{qpFunctionalCoherence}{lsCMatrix}(object, TFgenes, geneUniverse=rownames(object),
                                            chip, minRMsize=5, verbose=FALSE, clusterSize=1)
\S4method{qpFunctionalCoherence}{lspMatrix}(object, TFgenes, geneUniverse=rownames(object),
                                            chip, minRMsize=5, verbose=FALSE, clusterSize=1)
\S4method{qpFunctionalCoherence}{lsyMatrix}(object, TFgenes, geneUniverse=rownames(object),
                                            chip, minRMsize=5, verbose=FALSE, clusterSize=1)
\S4method{qpFunctionalCoherence}{matrix}(object, TFgenes, geneUniverse=rownames(object),
                                         chip, minRMsize=5, verbose=FALSE, clusterSize=1)
\S4method{qpFunctionalCoherence}{list}(object, geneUniverse=unique(c(names(object), unlist(object, use.names=FALSE))),
                                       chip, minRMsize=5, verbose=FALSE, clusterSize=1)
}
\arguments{
  \item{object}{object containing the transcriptional regulatory modules for which
       we want to estimate their functional coherence. It can be an adjacency matrix
       of the undirected graph representing the transcriptional regulatory network or
       a list of gene target sets where the name of the entry should be the transcription
       factor identifier.}
  \item{TFgenes}{when the input object is a matrix, it is required to provide a vector of
       transcription factor gene identifiers (which should match somewhere in the row and
       column names of the matrix.}
  \item{geneUniverse}{vector of all genes considered in the analysis. By default
       it equals the rows and column names of \code{object} when it is a matrix, or the
       set of all different gene identifiers occuring in \code{object} when it is a list.}
  \item{chip}{name of the \code{.db} package containing the Gene Ontology (GO)
       annotations.}
  \item{minRMsize}{minimum size of the target gene set in each regulatory module
       where functional enrichment will be calculated and thus where functional
       coherence will be estimated.}
  \item{verbose}{logical; if TRUE the function will show progress on the
       calculations; if FALSE the function will remain quiet (default).}
  \item{clusterSize}{size of the cluster of processors to employ if we wish to
       speed-up the calculations by performing them in parallel. A value of 1
       (default) implies a single-processor execution. The use of a cluster of
       processors requires having previously loaded the packages \code{snow}
       and \code{rlecuyer}.}
}
\details{
This function estimates the functional coherence of a transcriptional regulatory
network represented by means of an undirected graph encoded by an adjacency matrix
and of a set of transcription factor genes. The functional coherence of a
transcriptional regulatory network is calculated as specified by Castelo and
Roverato (2009) and corresponds to the distribution of individual functional
coherence values of every of the regulatory modules of the network each of them
defined as a transcription factor and its set of putatively regulated target
genes. In the calculation of the functional coherence value of a regulatory
module, Gene Ontology (GO) annotations are employed through the given annotation
\code{.db} package and the conditional hyper-geometric test implemented in the
\code{GOstats} package from Bioconductor.
}
\value{
A list with three slots, a first one containing the transcriptional regulatory
network as a list of regulatory modules and their targets, a second one
containing this same network but including only those modules with GO BP
annotations and a third one consisting of a vector of functional coherence values.
}
\references{
Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpAvgNrr}}
  \code{\link{qpGraph}}
}
\examples{
library(org.EcK12.eg.db)

# load RegulonDB data from this package
data(EcoliOxygen)

# pick two TFs from the RegulonDB data in this package

TFgenes <- c("mhpR", "iscR")

# get their Entrez Gene Identifiers
TFgenesEgIDs <- unlist(mget(TFgenes, AnnotationDbi::revmap(org.EcK12.egSYMBOL)))

# get all genes involved in their regulatory modules from
# the RegulonDB data in this package
mt <- match(filtered.regulon6.1[,"EgID_TF"], TFgenesEgIDs)

allGenes <- as.character(unique(as.vector(
            as.matrix(filtered.regulon6.1[!is.na(mt),
                                          c("EgID_TF","EgID_TG")]))))

mtTF <- match(filtered.regulon6.1[,"EgID_TF"],allGenes)
mtTG <- match(filtered.regulon6.1[,"EgID_TG"],allGenes)

# select the corresponding subset of the RegulonDB data in this package
subset.filtered.regulon6.1 <- filtered.regulon6.1[!is.na(mtTF) & !is.na(mtTG),]
TFi <- match(subset.filtered.regulon6.1[,"EgID_TF"], allGenes)
TGi <- match(subset.filtered.regulon6.1[,"EgID_TG"], allGenes)
subset.filtered.regulon6.1 <- cbind(subset.filtered.regulon6.1,
                                    idx_TF=TFi, idx_TG=TGi)

# build an adjacency matrix representing the transcriptional regulatory
# relationships from these regulatory modules
p <- length(allGenes)
adjacencyMatrix <- matrix(FALSE, nrow=p, ncol=p)
rownames(adjacencyMatrix) <- colnames(adjacencyMatrix) <- allGenes
idxTFTG <- as.matrix(subset.filtered.regulon6.1[,c("idx_TF","idx_TG")])
adjacencyMatrix[idxTFTG] <-
  adjacencyMatrix[cbind(idxTFTG[,2],idxTFTG[,1])] <- TRUE

# calculate functional coherence on these regulatory modules
fc <- qpFunctionalCoherence(adjacencyMatrix, TFgenes=TFgenesEgIDs,
                            chip="org.EcK12.eg.db")

print(sprintf("the \%s module has a FC value of \%.2f",
              mget(names(fc$functionalCoherenceValues),org.EcK12.egSYMBOL),
              fc$functionalCoherenceValues))
}
\keyword{models}
\keyword{multivariate}
