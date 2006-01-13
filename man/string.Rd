\name{fp.to.string}
\alias{fp.to.string}
\title{
  Generates a string representation of a binary fingerprint
}
\description{
    The function returns a string of 1's and 0's corresponding to the 
    fingerprint vector supplied        
}
\usage{
fp.to.string(fp, size=1024)
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
A string of 1's and 0's
}
\seealso{
    \code{\link{fp.from.bstring}}
}
\examples{
# make a fingerprint vector
fp <- fp.from.bstring("110011")

# print out the string representation
print(fp.to.string(fp, 6))
}
\keyword{logic}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
