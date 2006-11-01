\name{distance}
\alias{distance}
\title{
Calculates the Distance Between Two Fingerprints
}
\description{
    A number of distance metrics can be calculated for binary fingerprints.
    These metrics can be used to evaluate similarity/dissimilarity between
    fingerprints and hence are useful for clustering purposes. The function
    currently allows the evaluation of 4 distance metrics
    \itemize{
        \item Euclidean
        \item Tanimoto
        \item Dice
        \item Modified Tanimoto
    } The default metric is the Tanimoto coefficient.
    In the case of the last 3, the value is actually a similarity value and hence
    the distance metric is obtained by subtracting the obtained value from 1.0.    
}
\usage{
distance(fp1, fp2, method)
}
\arguments{
    \item{fp1}{
      An object of class \code{fingerprint}
    }
    \item{fp2}{
      An object of class \code{fingerprint}
    }
    \item{method}{
        The type of distance metric desired. Alternative values are
        \code{euclidean} and \code{dice} and \code{mt}. Partial matching
	is supported and the deault is \code{tanimoto}
    }
}
\value{
Numeric value representing the distance in the specified metric between the
supplied fingerprint objects
}
\examples{
# make a 2 fingerprint vectors
fp1 <- new("fingerprint", nbit=6, bits=c(1,2,5,6))
fp2 <- new("fingerprint", nbit=6, bits=c(1,2,5,6))

# calculate the tanimoto coefficient
distance(fp1,fp2) # should be 1

# Invert the second fingerprint
fp3 <- !fp2

distance(fp1,fp3) # should be 0
}
\references{Fligner, M.A.; Verducci, J.S.; Blower, P.E.;
A Modification of the Jaccard-Tanimoto Similarity Index for
Diverse Selection of Chemical Compounds Using Binary Strings,
\emph{Technometrics}, 2002, \emph{44}(2), 110-119
}
    
\keyword{logic}
\author{Rajarshi Guha \email{rguha@indiana.edu}}
