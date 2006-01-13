\name{fp.fold}
\alias{fp.fold}
\title{
  Fold a fingerprint
}
\description{
    In many situations a fingerprint is generated using a large length (such as 1024 bits or more).
    As a result of this, the fingerprints for a dataset can be very sparse. One approach to increasing 
    bit density of such fingerprints is to fold them. This is performed by dividing the original
    fingerprint bitstring into two substrings of equal length and then perform an XOR on
    the two substrings. 

    It should be noted that many fingerprint generating routines will perform this internally.
}
\usage{
fp.fold(fp, size=1024)
}
\arguments{
    \item{fp}{
        The fingerprint to fold. 
    }
    \item{size}{
        Length of the bit string represented by the fingerprint vector. Note that the size of
        should be an even number.    
    }
}
\value{
The folded fingerprint vector. The size of this vector should be taken as half that of
the original vector
}
\seealso{
    \code{\link{fp.xor}}}
}
\examples{
# make a fingerprint vector
fp <- fp.from.bstring("110011")
fp.fold( fp, size=6)
}
\keyword{logic}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
