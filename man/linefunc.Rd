\name{cdk.lf, moe.lf}
\alias{cdk.lf}
\alias{moe.lf}
\title{
    Functions to parse lines from fingerprint files
}
\description{
\code{cdk.lf} and \code{moe.lf} takes a single line and parses it to produce 
a vector of integers which represents the position of the 'on' bits in
a fingerprint. This allows the user to use \code{read.fp} with arbitrary fingerprint
files. A new file format can be handled by defining a new line parser function.
Currently the two functions process fingerprint files obtained from the 
CDK (\url{http://cdk.sourceforge.net}) and MOE (\url{http://chemcomp.com}).

}
\usage{
    cdk.lf(line)
    moe.lf(line)
}
\arguments{
    \item{line}{
        The line to parse
    }
}
\value{
A vector of integers representing 'on' bits 
}
\author{Rajarshi Guha (\email{rajarshi@presidency.com})}
\keyword{logic}



