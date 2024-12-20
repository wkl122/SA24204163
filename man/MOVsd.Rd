% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changepoint-detection.R
\name{MOVsd}
\alias{MOVsd}
\title{MOVsd}
\usage{
MOVsd(data, n)
}
\arguments{
\item{data}{the time series after diff}

\item{n}{the size of moving window}
}
\value{
the moving standard deviation
}
\description{
this function can compute the moving mean of time series
}
\examples{
\dontrun{
    x<-rbinom(20,size=10,prob=0.5)
    MOVsd(x,5)
}
}
