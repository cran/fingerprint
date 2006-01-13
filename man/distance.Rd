\name{fp.distance}
\alias{fp.distance}
\title{
  Calculates the distance between two fingerprint vectors
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
fp.distance(fp1, fp2, size=1024, type='tanimoto', ...)
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
    \item{type}{
        The type of distance metric desired. Alternative values are
        \code{euclidean} and \code{dice} and \code{mt}
    }
    \item{...}{
        Currently not used, but will be used to supply arguments
        to the Tversky metric (a generalization of the Tanimoto and Dice
       metrics)
    }
}
\value{
Numeric representing the distance in the specified metric between the
supplied fingerprint vectors
}
\examples{
# make a 2 fingerprint vectors
fp1 <- fp.from.bstring("110011")
fp2 <- fp.from.bstring("110011")

# calculate the tanimoto coefficient
fp.distance(fp1,fp2,6) # should be 1

# Invert the second fingerprint
fp3 <- fp.not(fp2, 6)

fp.distance(fp1,fp3,6) # should be 0
}
\references{Fligner, M.A.; Verducci, J.S.; Blower, P.E.;
A Modification of the Jaccard-Tanimoto Similarity Index for
Diverse Selection of Chemical Compounds Using Binary Strings,
\emph{Technometrics}, 2002, \emph{44}(2), 110-119
}
    
\keyword{logic}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
