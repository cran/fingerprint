\name{fp.and, fp.or, fp.not}
\alias{fp.and}
\alias{fp.or}
\alias{fp.not}
\alias{fp.xor}
\title{
    Logical operators for fingerprints
}
\description{
These functions perform logical operatiosn (AND, OR, NOT) on the supplied
binary fingerprints. Thus for two fingerprints A and B we have
\describe{
\item{fp.and}{Logical AND}
\item{fp.or}{Logical OR}
\item{fp.not}{Logical NOT (negation)}
\item{fp.xor}{Exclusive OR}
}
}
\usage{
fp.and(fp1, fp2, size=1024)
fp.or(fp1, fp2, size=1024)
fp.not(fp1, size=1024)
fp.xor(fp1, fp2, size=1024)
}
\arguments{
    \item{fp1}{
        A fingerprint vector
    }
    \item{fp2}{
        A fingerprint vector
    }
    \item{size}{
        The length of the fingerprints being considered
    }
}
\value{
A fingerprint vector
}
\examples{
# make 2 fingerprints
fp1 <- fp.from.bstring("11001011")
fp2 <- fp.from.bstring("10011010")
and <- fp.and(fp1,fp2, size=8)
or <- fp.or(fp1,fp2, size=8)
not <- fp.not(fp1, size=8)

fp.to.string(and,size=8)
fp.to.string(or,size=8)
fp.to.string(not,size=8)
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}



