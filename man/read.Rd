\name{fp.read, fp.read.to.matrix}
\alias{fp.read}
\alias{fp.read.to.matrix}
\title{
    Functions to read fingerprints from files
}
\description{
\code{fp.read} reads in a set of fingerprints from a file. Fingerprint 
output from the CDK, MOE and BCI can be handled.

Each fingerprint is represented as an integer vector whose elements indicate which
positions of the bit string are set to 1. Since the fingerprint is simply a \code{numeric()}
most functions will require the length of the actual bit string to be passed in.
\code{fp.read} returns a \code{list} structure, each element being a fingerprint vector.

\code{fp.read.to.matrix} is a utility function that reads the fingerprints directly to
matrix form (columns are the bit positions and the rows are the objects whose fingerprints
have been evaluated)
}
}
\usage{
fp.read(f='fingerprint.txt', lf=cdk.lf, header=FALSE)
fp.read.to.matrix(f='fingerprint.txt', size=1024, lf=cdk.lf, header=FALSE)
}
\arguments{
  \item{f}{
    File containing the fingperprints
  }
  \item{size}{
    The length of the fingerprints being considered
  }
  \item{lf}{
    A line reading function that parses a single line from
    a fingerprint file. Currently, two such functions are provided
    that parse the fingerprints from the output of the CDK and MOE
    respectively. See cdk.lf and moe.lf
  }
  \item{header}{
    Indicates whether the first line of the fingerprint file is
    a header line
  }
}
\seealso{
  \code{\link{cdk.lf}},
  \code{\link{moe.lf}},
  \code{\link{bci.lf}}
}
\value{
  A \code{list} or \code{matrix} of fingerprint vectors
}
\author{Rajarshi Guha \email{rguha@indiana.edu}}
\keyword{logic}



