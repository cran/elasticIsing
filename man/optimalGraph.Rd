\name{optimalGraph}
\alias{optimalGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal graph selection
}
\description{
Selects the optimal graph from the results of \code{\link{elasticIsing}}
}
\usage{
optimalGraph(object, cost)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An \code{elasticIsing} object, resulting from  \code{\link{elasticIsing}}
}
  \item{cost}{
Cost function to use (from the cvTools package). One of \code{'mspe'}, \code{'rmspe'}, \code{'mape'}, \code{'tmspe'}. or \code{'rtmspe'}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
