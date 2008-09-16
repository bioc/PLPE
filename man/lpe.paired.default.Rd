\name{lpe.paired.default}
\alias{lpe.paired.default}
\title{Local Pooled Error Test for Paired Data}
\description{
This invetigates differential expression for paired high-throughput data.
}

\usage{
\method{lpe.paired}{default}(x, design, data.type, q=0.01, probe.ID = NULL, estimator="median", w=0.5, w.estimator="fixed", iseed=1234, ...)  
}

\arguments{
   \item{x}{data matrix}
   \item{design}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{q}{quantile for intervals of intensities}
   \item{probe.ID}{probe set IDs; if NULL, row numbers are assigned.}
   \item{data.type}{data type: 'ms' for mass spectrometry data, 'cdna' for cDNA microarray data }
   \item{estimator}{specification for the estimator: 'median', 'mean' and 'huber' }
   \item{w}{weight paramter between individual variance estimate and pooling variance estimate, 0<= w <=1}
   \item{w.estimator}{two approaches to estimate the weight: 'random' or 'fixed' }
   \item{iseed}{seed number}
   \item{...}{other arguments}
}

\value{
   \item{design}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{data.type}{data type: 'ms' for mass spectrometry data, 'cdna' for cDNA microarray data }
   \item{q}{quantile for intervals of intensities}
   \item{estimator}{specification for the estimator: 'median', 'mean' and 'huber'}
   \item{w.estimator}{two approaches to estimate the weight: 'random' or 'fixed'}
   \item{w}{weight paramter between individual variance estimate and pooling variance estimate, 0<= w <=1}
   \item{test.out}{matrix for test results}
}

\references{
Cho H, Smalley DM, Ross MM, Theodorescu D, Ley K and Lee JK (2007).  Statistical Identification of Differentially Labelled Peptides from Liquid Chromatography Tandem Mass Spectrometry, Proteomics, 7:3681-3692.
 }


\author{
HyungJun Cho and Jae K. Lee 
}

\seealso{
    \code{\link{lpe.paired}}
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
summary(out)
}

\keyword{models}


