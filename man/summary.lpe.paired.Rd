\name{summary.lpe.paired}
\alias{summary.lpe.paired}
\title{Local Pooled Error Test for Paired Data}
\description{
This summarize lpe.paired.
}
\usage{
\method{summary}{lpe.paired}(object,...)
}

\arguments{
   \item{x}{an object created from  lpe.paired}
   \item{...} {other arguments}
}


\references{
Cho H, Smalley DM, Ross MM, Theodorescu D, Ley K and Lee JK (2007).  Statistical Identification of Differentially Labelled Peptides from Liquid Chromatography Tandem Mass Spectrometry, Proteomics, 7:3681-3692.
 }


\author{
HyungJun Cho and Jae K. Lee 
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
summary(out)
}

\keyword{internal}


