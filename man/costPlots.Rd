\name{costPlots}
\alias{costPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D plots of cost functions
}
\description{
This function creates a PDF file with the 3D plots of cost functions given alpha and lambda. By default, accuracy is shown, which is defined as -cost.
}
\usage{
costPlots(object, filename = "elasticIsing.pdf", width = 8, 
    height = 5, theta = 25, phi = 30, ticktype = "simple", 
    accuracy = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An \code{elasticIsing} object, resulting from  \code{\link{elasticIsing}}
}
  \item{filename}{
Name of the PDF file to create
}
  \item{width}{
Width of the PDF file to creare
}
  \item{height}{
Height of the PDF file to create
}
  \item{theta}{
See \code{\link{persp}}
}
  \item{phi}{
See \code{\link{persp}}
}
  \item{ticktype}{
See \code{\link{persp}}
}
  \item{accuracy}{
See \code{\link{persp}}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
