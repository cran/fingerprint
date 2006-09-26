\name{fp.num}
\alias{fp.num}
\title{
    The number of bits that are on in a fingerprint
}
\description{
}
\usage{
fp.num(fp)
}
\arguments{
\item{fp}{
    A fingerprint vector
}
}
\value{
An integer indicating how many bits are on
}
\examples{
fp <- fp.from.bstring("11001011")
fp.num(fp)
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}


