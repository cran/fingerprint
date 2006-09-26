\name{fp.from.bstring}
\alias{fp.from.bstring}
\title{
  Generates a fingerprint vector from a string representation of the binary
  fingerprint
}
\description{
    The function accepts a string containing 0's or 1's and generates a 
    fingerprint vector that can be used with the other functions in this package.
    An important feature is that the resultant vector has no representation of
    the length of the fingerprint
}
\usage{
fp.from.bstring(fp)
}
\arguments{
\item{fp}{
    A string of 1's and 0's with no spaces
}
}
\value{
A vector of integers representing positions that contain a 1
}
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}
