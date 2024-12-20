% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changepoint-detection.R
\name{changepoint}
\alias{changepoint}
\title{changepoint}
\usage{
changepoint(D, ni, increase, m = 250)
}
\arguments{
\item{D}{the time series}

\item{ni}{the size of moving window}

\item{increase}{detect increase changepoints with increase=1 or decrease changepoints with increase=2}

\item{m}{the size of window}
}
\value{
a list containS BDM,movemean,movesd,BDI,BDIquant and test of the time series
}
\description{
this function can detect changepoints of time series using sliding visibility graph algorithm and backward degree metric
}
\examples{
\dontrun{
         set.seed(12)
         x<-c()
         for(i in 1:90){
           x<-c(x,rnorm(1,0,0.3))
         }
         for(i in 91:99){
           x<-c(x,((i-90)/6)^(4)+rnorm(1,0,0.3))
         }
         for(i in 100:111){
           x<-c(x,((112-i)/8)^4+rnorm(1,0,0.3))
         }
         for(i in 112:200){
           x<-c(x,rnorm(1,0,0.3))
         }
         changepoint(x,9,1)
         changepoint(x,12,2)
}
}
