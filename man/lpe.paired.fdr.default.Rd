\name{lpe.paired.fdr.default}
\alias{lpe.paired.fdr.default}
\title{FDR for PLPE}
\description{
This computes FDR for PLPE.
}

\usage{
\method{lpe.paired.fdr}{default}(x, obj, n.iter=5, lambda=0.9, ...) 
}

\arguments{
   \item{x}{data matrix}
   \item{obj}{object created from  lpe.paired}
   \item{n.iter}{number of iterations}
   \item{lambda}{numeric vector of probabilities with values in [0,1]}
   \item{...}{other argument}
}

\value{
   \item{design}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{data.type}{data type: 'ms' for mass spectrometry data, 'cdna' for cDNA microarray data }
   \item{estimator}{specification for the estimator: 'median', 'mean' and 'huber'}
   \item{w.estimator}{two approaches to estimate the weight: 'random' or 'fixed'}
   \item{w}{weight paramter between individual variance estimate and pooling variance estimate, 0<= w <=1}
   \item{pi0}{estimated proportion of non-null peptides}
   \item{FDR}{matrix for test results including FDRs}
   \item{...}{other arguments}
}

\references{
Cho H, Smalley DM, Ross MM, Theodorescu D, Ley K and Lee JK (2007).  Statistical Identification of Differentially Labelled Peptides from Liquid Chromatography Tandem Mass Spectrometry, Proteomics, 7:3681-3692.
 }


\author{
HyungJun Cho and Jae K. Lee 
}

\seealso{
    \code{\link{lpe.paired.fdr}}
}

\examples{

#LC-MS/MS proteomic data for platelets MPs
library(PLPE)
data(plateletSet)
x <- exprs(plateletSet)
x <- log2(x) 

cond <- c(1, 2, 1, 2, 1, 2)
pair <- c(1, 1, 2, 2, 3, 3)
design <- cbind(cond, pair)

out <- lpe.paired(x, design, q=0.1, data.type="ms")
out.fdr <- lpe.paired.fdr(x,obj=out)
out.fdr$FDR[1:10,]

}

\keyword{models}


