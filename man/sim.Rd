\name{fp.sim.matrix}
\alias{fp.sim.matrix}
\title{
  Calculates the similarity matrix for a set of fingerprints
}
\description{
Given a set of fingerprints, a pairwise similarity can be calculated using the
various distance metrics defined for binary strings. This function calculates
the pairwise similarity matrix for a set of fingerprint vectors supplied in a \code{list} 
structure. Any of the distance metrics provided by \code{fp.distance} can be used and the 
default is the Tanimoto metric.

Note that if the the Euclidean distance is specified then the resultant matrix is a
distance matrix and not a similarity matrix
}
\usage{
fp.sim.matrix(fplist, size=1024, type='tanimoto')
}
\arguments{
    \item{fplist}{
        A list structure with each element being a fingerprint vector. This
        structure will generally be obtained from \code{fp.read}.
    }
    \item{size}{
        Length of the bit string represented by the fingerprint vector
    }
    \item{type}{
    The type of distance metric to use. Alternatives are \code{euclidean} and
    \code{dice} and \code{mt}
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

fp.sim.matrix( list(fp1,fp2,fp3), size=6)
}
\keyword{logic}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
