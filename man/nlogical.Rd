\name{fpNOr, fpNAnd, fpNOff}
\alias{fpNOr}
\alias{fpNAnd}
\alias{fpNOff}
\title{
    These functions evaluate the number of bits in various conditions
}
\description{
These functions evaluate the number of bits common (i.e., set or unset)
between two fingerprints. These are mainly used to evaluate the various 
distance metrics. So for two fingerprints A and B we have
\describe{
\item{fpNOr}{Number of bits on in A but not on in B}
\item{fpNAnd}{Number of bits on in A that are also on in B}
\item{fpNOff}{Number of bits off in A that are also off in B}
}
}
\usage{
fpNOr(fp1, fp2, size=1024)
fpNAnd(fp1, fp2, size=1024)
fpNOff(fp1, fp2, size=1024)
}
\arguments{
\item{fp1}{
    A fingerprint vector
}
\item{fp2}{
    A fingerprint vector
}
\item{size}{
    The length of the fingerprints represented by the fingerprint vectors
}
}
\value{
    An integer indicating the number of bits for the specific condition
}
\examples{
# make 2 fingerprints
fp1 <- fp.from.bstring("11001011")
fp2 <- fp.from.bstring("11010100")
fpNOr(fp1,fp2,8)
fpNAnd(fp1,fp2,8)
fpNOff(fp1,fp2,8)
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}


