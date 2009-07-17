\name{qpFunctionalCoherence}
\alias{qpFunctionalCoherence}

\title{
Functional coherence estimation
}
\description{
Estimates functional coherence for a given transcriptional regulatory network.
}
\usage{
qpFunctionalCoherence(I, TFgenes, chip, minRMsize=5, verbose=FALSE)
}
\arguments{
  \item{I}{incidence matrix of the undirected graph representing the
       transcriptional regulatory network.}
  \item{TFgenes}{vector of transcription factor gene (matching the genes at the
       rows and column names of \code{I}.}
  \item{chip}{name of the \code{.db} package containing the Gene Ontology (GO)
       annotations.}
  \item{minRMsize}{minimum size of the target gene set in each regulatory module
       where functional enrichment will be calculated and thus where functional
       coherence will be estimated.}
  \item{verbose}{logical; if TRUE the function will show progress on the
       calculations; if FALSE the function will remain quiet (default).}
}
\details{
This function estimates the functional coherence of a transcriptional regulatory
network represented by means of an undirected graph encoded by an incidence matrix
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
library(annotate)
library(org.EcK12.eg.db)

# load RegulonDB data from this package
data(EcoliOxygen)

# pick two TFs from the RegulonDB data in this package

TFgenes <- c("mhpR", "iscR")

# get their Entrez Gene Identifiers
TFgenesEgIDs <- unlist(mget(TFgenes, revmap(getAnnMap("SYMBOL", "org.EcK12.eg.db"))))

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

# build an incidence matrix representing the transcriptional regulatory
# relationships from these regulatory modules
p <- length(allGenes)
incidenceMatrix <- matrix(FALSE, nrow=p, ncol=p)
rownames(incidenceMatrix) <- colnames(incidenceMatrix) <- allGenes
idxTFTG <- as.matrix(subset.filtered.regulon6.1[,c("idx_TF","idx_TG")])
incidenceMatrix[idxTFTG] <-
  incidenceMatrix[cbind(idxTFTG[,2],idxTFTG[,1])] <- TRUE

# calculate functional coherence on these regulatory modules
fc <- qpFunctionalCoherence(incidenceMatrix, TFgenesEgIDs, "org.EcK12.eg.db")

print(sprintf("the \%s module has a FC value of \%.2f",
              mget(names(fc$functionalCoherenceValues),org.EcK12.egSYMBOL),
              fc$functionalCoherenceValues))
}
\keyword{models}
\keyword{multivariate}
