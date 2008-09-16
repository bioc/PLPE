\name{lpe.paired.fdr}
\alias{lpe.paired.fdr}
\title{FDR for PLPE}
\description{
This computes FDR for PLPE.
}

\usage{
lpe.paired.fdr(x,...)
}

\arguments{
   \item{x}{data matrix}
   \item{...}{other arguments}
}

\references{
Cho H, Smalley DM, Ross MM, Theodorescu D, Ley K and Lee JK (2007).  Statistical Identification of Differentially Labelled Peptides from Liquid Chromatography Tandem Mass Spectrometry, Proteomics, 7:3681-3692.
 }


\author{
HyungJun Cho and Jae K. Lee 
}

\seealso{
    \code{\link{lpe.paired.fdr.default}}
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


