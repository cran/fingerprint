\name{fp.to.vector}
\alias{fp.to.vector}
\title{
  Converts the fingerprint from a nD cube representation to an
  nD sphere representation
}
\description{
 Ordinarily, a binary fingerprint can be considered to represent a 
 corner of a nD hypercube. However in many cases using such a representation
 can lead to a very sparse space. Consequently one approach is to convert
 the fingerprint so that it represents points on a nD unit hypersphere.

 The resultant fingerprint is then a nD coordinate.
}
\usage{
fp.to.vector(fp, size=1024)
}
\arguments{
    \item{fp}{
        A fingerprint vector
    }
    \item{size}{
        Length of the bit string represented by the fingerprint vector
    }
}
\value{
A numeric of length equal to size corresponding to a unit vector for a point
on the nD hypersphere
}
\examples{
# make a fingerprint vector
fp <- fp.from.bstring("110011")

# print out the string representation
vec <- fp.to.vector(fp,6)
}
\keyword{logic}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
