\name{generate.null}
\alias{generate.null}
\title{Local Pooled Error Test for Paired Data}
\description{
This is internal function in lpe.paired.
}

\usage{
generate.null(x, design, q)
}

\arguments{
   \item{x}{data matrix}
   \item{design}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{q}{quantile for intervals of intensities}
}

\value{
   \item{design}{design matrix; condition index in the first column and pair index in the sceond column}
   \item{q}{quantile for intervals of intensities}
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

\keyword{internal}


