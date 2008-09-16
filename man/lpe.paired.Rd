\name{lpe.paired}
\alias{lpe.paired}
\title{Local Pooled Error Test for Paired Data}
\description{
This invetigates differential expression for paired high-throughput data.
}
\usage{
lpe.paired(x,...)
}

\arguments{
   \item{x}{an object for which the extraction of model lpe.paired is meaningful.}
   \item{...}{other arguments}
}

\value{
   \item{x}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{...}{data type: 'ms' for mass spectrometry data, 'cdna' for cDNA microarray data }
}

\references{
Cho H, Smalley DM, Ross MM, Theodorescu D, Ley K and Lee JK (2007).  Statistical Identification of Differentially Labelled Peptides from Liquid Chromatography Tandem Mass Spectrometry, Proteomics, 7:3681-3692.
 }


\author{
HyungJun Cho and Jae K. Lee 
}

\seealso{
    \code{\link{lpe.paired.default}}
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
out$test.out[1:10,]

}

\keyword{models}


