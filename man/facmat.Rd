\name{fp.factor.matrix}
\alias{fp.factor.matrix}
\title{
  Converts a list of fingerprints to a data.frame of factors
}
\description{
This function will convert  a \code{list} of fingerprint
to a \code{data.frame} of factors with levels 1 and 0.
}
\usage{
fp.factor.matrix(fplist, size=1024)
}
\arguments{
    \item{fplist}{
        A list structure with each element being a fingerprint vector. This
        structure will generally be obtained from \code{fp.read}.
    }
    \item{size}{
        Length of the bit string represented by the fingerprint vector
    }
}
\value{
A matrix with dimensions equal to \code{(length(fplist), length(fplist))}
}
\seealso{
    \code{\link{fp.distance}}, \code{\link{fp.read}}
}
\examples{
# make a fingerprint vector
fp1 <- fp.from.bstring("110011")
fp2 <- fp.from.bstring("100111")
fp3 <- fp.from.bstring("011111")

fp.factor.matrix( list(fp1,fp2,fp3), size=6)
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}
