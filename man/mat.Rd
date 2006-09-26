\name{fp.to.matrix}
\alias{fp.to.matrix}
\title{
  Converts a list structure of fingerpring vectors to a matrix
}
\description{
    In general, fingerprint data is read from a file or obtained via
    calls to an external generator and the return value is a list of fingerprints.
    This function takes the list and returns a matrix having number of rows equal to 
    the number of fingerprints and the number of columns equal to the length of
    the fingerprint. Each element is 1 or 0 (1's being specified by the positions
    in each fingerprint vector)
}
\usage{
fp.to.matrix(fplist, size=1024)
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
A matrix with dimensions equal to \code{(length(fplist), size)}
}
\seealso{
    \code{\link{fp.distance}}, \code{\link{fp.read}}
}
\examples{
# make a fingerprint vector
fp1 <- fp.from.bstring("110011")
fp2 <- fp.from.bstring("100111")
fp3 <- fp.from.bstring("011111")

fp.to.matrix( list(fp1,fp2,fp3), size=6)
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}
