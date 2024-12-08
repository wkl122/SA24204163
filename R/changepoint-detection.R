library(zoo)
library(ggplot2)
library(cpm)
library(Rcpp)
#' @title visibility
#' @description this function can transform time series into a complex network using sliding visibility graph algorithm
#' @param data the time series after diff
#' @param m the size of the window
#' @return a graph
#' @examples
#' \dontrun{
#'     x<-rnorm(100,0,1)
#'     visibility(x,100)
#' }
#' @export
visibility<-function(data,m){
  graph<-matrix(0,length(data),length(data))
  if(length((data))==2){
    graph<-matrix(c(0,1,1,0),2,2)
    return(graph)
  }
  for (i in 1:(length(data)-2)) {
    for(j in (i+2):(min((i+2+m),length(data)))){
      a<-0
      for(k in (i+1):(j-1)){
        b<-(data[i]+(data[j]-data[i])*(k-i)/(j-i))
        if((data[k]<b)){
          a<-a
        }else{
          a<-a+1
        }
      }
      if(a==0){
        graph[i,j]<-1
      }else{
        graph[i,j]<-0
      }
    }

  }
  for(i in 1:(length(data)-1)){
    graph[i,i+1]<-1
  }
  graph<-graph+t(graph)
  return(graph)
}
sourceCpp("..//SA24204163//src//inst//BDM.cpp")
# BDM<-function(graph){
#   c<-rep(0,length(graph[1,]))
#   for(i in 1:length(graph[1,])){
#     for(j in 1:i){
#       c[i]<-c[i]+graph[i,j]
#     }
#   }
#   return(c)
# }

#' @title MOVmean
#' @description this function can compute the moving mean of time series
#' @param data the time series after diff
#' @param n the size of moving window
#' @return the moving mean
#' @examples
#' \dontrun{
#'     x<-rbinom(20,size=10,prob=0.5)
#'     MOVmean(x,5)
#' }
#' @export
MOVmean<-function(data,n){
  mean<-rep(0,length(data))
  if(n==1){
    mean<-data
  }else if(n==2){
    for(i in 1:(n/2)){
      mean[i]<-data[i]
    }
    for(i in (n/2+1):(length(data)-n/2+1)){
      mean[i]<-mean(data[(i-n/2):(i+n/2-1)])
    }
  }else if((n>2)&(n%%2==0)){
    for(i in 1:(n/2)){
      mean[i]<-mean(data[1:(i+n/2-1)])
    }
    for(i in (n/2+1):(length(data)-n/2+1)){
      mean[i]<-mean(data[(i-n/2):(i+n/2-1)])
    }
    for(i in (length(data)-n/2+2):length(data)){
      mean[i]<-mean(data[(i-n/2):length(data)])
    }
  }else {
    n2<-floor(n/2)
    n1<-ceiling(n/2)
    for(i in 1:n2){
      mean[i]<-mean(data[1:(i+n2)])
    }
    for(i in (n2+1):(length(data)-n2)){
      mean[i]<-mean(data[(i-n2):(i+n2)])
    }
    for(i in (length(data)-n2+1):length(data)){
      mean[i]<-mean(data[(i-n2):length(data)])
    }
  }
  return(mean)
}

#' @title MOVsd
#' @description this function can compute the moving mean of time series
#' @param data the time series after diff
#' @param n the size of moving window
#' @return the moving standard deviation
#' @examples
#' \dontrun{
#'     x<-rbinom(20,size=10,prob=0.5)
#'     MOVsd(x,5)
#' }
#' @export
MOVsd<-function(data,n){
  sd<-rep(0,length(data))
  if(n==1){
    sd<-rep(0,length(data))
  }else if(n==2){
    for(i in 1:(n/2)){
      sd[i]<-0
    }
    for(i in (n/2+1):(length(data)-n/2+1)){
      sd[i]<-sd(data[(i-n/2):(i+n/2-1)])
    }
  }else if((n>2)&(n%%2==0)){
    for(i in 1:(n/2)){
      sd[i]<-sd(data[1:(i+n/2-1)])
    }
    for(i in (n/2+1):(length(data)-n/2+1)){
      sd[i]<-sd(data[(i-n/2):(i+n/2-1)])
    }
    for(i in (length(data)-n/2+2):length(data)){
      sd[i]<-sd(data[(i-n/2):length(data)])
    }
  }else{
    n2<-floor(n/2)
    n1<-ceiling(n/2)
    for(i in 1:n2){
      sd[i]<-sd(data[1:(i+n2)])
    }
    for(i in (n2+1):(length(data)-n2)){
      sd[i]<-sd(data[(i-n2):(i+n2)])
    }
    for(i in (length(data)-n2+1):length(data)){
      sd[i]<-sd(data[(i-n2):length(data)])
    }
  }
  return(sd)
}

#' @title changepoint
#' @description this function can detect changepoints of time series using sliding visibility graph algorithm and backward degree metric
#' @param D the time series
#' @param ni the size of moving window
#' @param increase detect increase changepoints with increase=1 or decrease changepoints with increase=2
#' @param m the size of window
#' @return a list containS BDM,movemean,movesd,BDI,BDIquant and test of the time series
#' @examples
#' \dontrun{
#'          set.seed(12)
#'          x<-c()
#'          for(i in 1:90){
#'            x<-c(x,rnorm(1,0,0.3))
#'          }
#'          for(i in 91:99){
#'            x<-c(x,((i-90)/6)^(4)+rnorm(1,0,0.3))
#'          }
#'          for(i in 100:111){
#'            x<-c(x,((112-i)/8)^4+rnorm(1,0,0.3))
#'          }
#'          for(i in 112:200){
#'            x<-c(x,rnorm(1,0,0.3))
#'          }
#'          changepoint(x,9,1)
#'          changepoint(x,12,2)
#' }
#' @export
changepoint<-function(D,ni,increase,m=250){
  z<-D
  D<-diff(D)
  if(increase==1){
    for(i in 1:length(D)){
      if(D[i]<0){
        D[i]=-D[i]
      }
    }
  }else{
    for(i in 1:length(D)){
      if(D[i]>0){
        D[i]=0
      }else{
        D[i]=-D[i]
      }
    }
  }
  #plot(D,type='l')
  Net_test<-visibility (D,m)
  #print(Net_test)
  #backward network degree
  Net_Degree_test=BDM(Net_test)
  #plot(Net_Degree_test,type = 'l')
  #moving mean of backward degree time series
  MA2test=MOVmean(Net_Degree_test,ni)
  #moving std of backward degree time series
  MStd2test=MOVsd(Net_Degree_test,ni)
  #BDI
  MA2MStd2test=MA2test+MStd2test
  BDIquant<-(mean(MA2MStd2test)+2*sd(MA2MStd2test))
  t<-rep(0,length(MA2MStd2test))
  for(i in 1:length(MA2MStd2test)){
    if(MA2MStd2test[i]>BDIquant){
      t[i]<-1
    }
  }
  data1<-data.frame("time1"=1:length(z),z,"test1"=c(t,0))
  p<-ggplot(data=data1,aes(x=data1$time1, y = z,color=data1$test1)) + geom_point() +
    geom_line()+theme_bw()+ggtitle("changepoint") + theme(plot.title = element_text(hjust = 0.5))
  list<-list(Net_Degree_test,MA2test,MStd2test,MA2MStd2test,BDIquant,t,p)
  names(list)<-c("BDM","movmean","movsd","BDI","BDIquant","test","figure")
  return(list)
}
#changepoint(DSPtest,16,1)
#changepoint(DSPtest,23,2)
# set.seed(12)
# x<-c()
# for(i in 1:90){
#   x<-c(x,rnorm(1,0,0.3))
# }
# for(i in 91:99){
#   x<-c(x,((i-90)/6)^(4)+rnorm(1,0,0.3))
# }
# for(i in 100:111){
#   x<-c(x,((112-i)/8)^4+rnorm(1,0,0.3))
# }
# for(i in 112:200){
#   x<-c(x,rnorm(1,0,0.3))
# }
# changepoint(x,9,1)
#' @import Rcpp
#' @import microbenchmark
#' @import DAAG
#' @import twosamples
#' @import boot
#' @import MASS
#' @import mvtnorm
#' @import lpSolve
#' @import bootstrap
#' @import zoo
#' @import ggplot2
#' @importFrom stats sd
#' @import cpm
#' @import igraph
#' @useDynLib SA24204163
NULL
