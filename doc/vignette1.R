## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=F-------------------------------------------------------------------
library(ggplot2)
library(Rcpp)
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
sourceCpp("../src/inst/BDM.cpp")
# BDM<-function(graph){
#   c<-rep(0,length(graph[1,]))
#   for(i in 1:length(graph[1,])){
#     for(j in 1:i){
#       c[i]<-c[i]+graph[i,j]
#     }
#   }
#   return(c)
# }


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
  p<-ggplot(data=data1,aes(x=time1, y = z,color=test1)) + geom_point() +
    geom_line()+theme_bw()+ggtitle("changepoint") + theme(plot.title = element_text(hjust = 0.5))
  list<-list(Net_Degree_test,MA2test,MStd2test,MA2MStd2test,BDIquant,t,p)
  names(list)<-c("BDM","movmean","movsd","BDI","BDIquant","test","figure")
  return(list)
}

## ----eval=FALSE---------------------------------------------------------------
#  visibility<-function(data,m){
#    graph<-matrix(0,length(data),length(data))
#    if(length((data))==2){
#      graph<-matrix(c(0,1,1,0),2,2)
#      return(graph)
#    }
#    for (i in 1:(length(data)-2)) {
#      for(j in (i+2):(min((i+2+m),length(data)))){
#        a<-0
#        for(k in (i+1):(j-1)){
#          b<-(data[i]+(data[j]-data[i])*(k-i)/(j-i))
#          if((data[k]<b)){
#            a<-a
#          }else{
#            a<-a+1
#          }
#        }
#        if(a==0){
#          graph[i,j]<-1
#        }else{
#          graph[i,j]<-0
#        }
#      }
#  
#    }
#    for(i in 1:(length(data)-1)){
#      graph[i,i+1]<-1
#    }
#    graph<-graph+t(graph)
#    return(graph)
#  }

## ----fig.width=6, fig.height=4------------------------------------------------
library(ggplot2)
x<-c(0.87,0.49,0.36,0.83,0.87,0.49,0.36,0.83,0.87,0.49,0.36,0.83,0.87,0.49,0.36,0.83,0.87,0.49,0.36,0.83)
data <- data.frame(id = 1:20, x =x)
p <- ggplot(data, aes(x = id, y = x)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal()
segments1 <- data.frame(
  x = data$id[1:19],
  y = data$x[1:19], 
  xend = data$id[2:20],
  yend = data$x[2:20]  
)
segments2<-data.frame(
  x=data$id[c(1,1,1,2,5,5,5,6,9,9,9,10,13,13,13,14,17,17,18)],
  y=data$x[c(1,1,1,2,5,5,5,6,9,9,9,10,13,13,13,14,17,17,18)],
  xend=data$id[c(3,4,5,4,7,8,9,8,11,12,13,12,15,16,17,16,19,20,20)],
  yend=data$x[c(3,4,5,4,7,8,9,8,11,12,13,12,15,16,17,16,19,20,20)])
p <- p + geom_segment(data = segments1, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.8, color = "#CC0033")
p <- p + geom_segment(data = segments2, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.8, color = "#CC0033")
p <- p + theme_bw() + xlab("时间")
print(p)


## ----fig.width=6, fig.height=4------------------------------------------------
library(igraph)
g<-graph_from_adjacency_matrix(visibility(x,100), mode = "undirected")
V(g)$color <- c(rep(c("#FF66FF","#88c4e8","#FF6699","#33CC66"),5))  
E(g)$color <- "steelblue"                 
plot(g, vertex.label = V(g)$name, vertex.size =20, vertex.label.cex =1)

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('NumericVector BDM(NumericMatrix graph) {
   int m = graph.ncol();
   NumericVector backward_degree(m, 0.0);  

   for (int i = 0; i < m; ++i) {
     for (int j = 0; j <= i; ++j) { 
       backward_degree[i] += graph(i, j);
     }
   }
   return backward_degree;
}')
print(BDM(visibility(x,100)))

## ----eval=FALSE---------------------------------------------------------------
#  MOVmean<-function(data,n){
#    mean<-rep(0,length(data))
#    if(n==1){
#      mean<-data
#    }else if(n==2){
#      for(i in 1:(n/2)){
#        mean[i]<-data[i]
#      }
#      for(i in (n/2+1):(length(data)-n/2+1)){
#        mean[i]<-mean(data[(i-n/2):(i+n/2-1)])
#      }
#    }else if((n>2)&(n%%2==0)){
#      for(i in 1:(n/2)){
#        mean[i]<-mean(data[1:(i+n/2-1)])
#      }
#      for(i in (n/2+1):(length(data)-n/2+1)){
#        mean[i]<-mean(data[(i-n/2):(i+n/2-1)])
#      }
#      for(i in (length(data)-n/2+2):length(data)){
#        mean[i]<-mean(data[(i-n/2):length(data)])
#      }
#    }else {
#      n2<-floor(n/2)
#      n1<-ceiling(n/2)
#      for(i in 1:n2){
#        mean[i]<-mean(data[1:(i+n2)])
#      }
#      for(i in (n2+1):(length(data)-n2)){
#        mean[i]<-mean(data[(i-n2):(i+n2)])
#      }
#      for(i in (length(data)-n2+1):length(data)){
#        mean[i]<-mean(data[(i-n2):length(data)])
#      }
#    }
#    return(mean)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  MOVsd<-function(data,n){
#    sd<-rep(0,length(data))
#    if(n==1){
#      sd<-rep(0,length(data))
#    }else if(n==2){
#      for(i in 1:(n/2)){
#        sd[i]<-0
#      }
#      for(i in (n/2+1):(length(data)-n/2+1)){
#        sd[i]<-sd(data[(i-n/2):(i+n/2-1)])
#      }
#    }else if((n>2)&(n%%2==0)){
#      for(i in 1:(n/2)){
#        sd[i]<-sd(data[1:(i+n/2-1)])
#      }
#      for(i in (n/2+1):(length(data)-n/2+1)){
#        sd[i]<-sd(data[(i-n/2):(i+n/2-1)])
#      }
#      for(i in (length(data)-n/2+2):length(data)){
#        sd[i]<-sd(data[(i-n/2):length(data)])
#      }
#    }else{
#      n2<-floor(n/2)
#      n1<-ceiling(n/2)
#      for(i in 1:n2){
#        sd[i]<-sd(data[1:(i+n2)])
#      }
#      for(i in (n2+1):(length(data)-n2)){
#        sd[i]<-sd(data[(i-n2):(i+n2)])
#      }
#      for(i in (length(data)-n2+1):length(data)){
#        sd[i]<-sd(data[(i-n2):length(data)])
#      }
#    }
#    return(sd)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  changepoint<-function(D,ni,increase,m=250){
#    z<-D
#    D<-diff(D)
#    if(increase==1){
#      for(i in 1:length(D)){
#        if(D[i]<0){
#          D[i]=-D[i]
#        }
#      }
#    }else{
#      for(i in 1:length(D)){
#        if(D[i]>0){
#          D[i]=0
#        }else{
#          D[i]=-D[i]
#        }
#      }
#    }
#    #plot(D,type='l')
#    Net_test<-visibility (D,m)
#    #print(Net_test)
#    #backward network degree
#    Net_Degree_test=BDM(Net_test)
#    #plot(Net_Degree_test,type = 'l')
#    #moving mean of backward degree time series
#    MA2test=MOVmean(Net_Degree_test,ni)
#    #moving std of backward degree time series
#    MStd2test=MOVsd(Net_Degree_test,ni)
#    #BDI
#    MA2MStd2test=MA2test+MStd2test
#    BDIquant<-(mean(MA2MStd2test)+2*sd(MA2MStd2test))
#    t<-rep(0,length(MA2MStd2test))
#    for(i in 1:length(MA2MStd2test)){
#      if(MA2MStd2test[i]>BDIquant){
#        t[i]<-1
#      }
#    }
#    data1<-data.frame("time1"=1:length(z),z,"test1"=c(t,0))
#    p<-ggplot(data=data1,aes(x=time1, y = z,color=test1)) + geom_point() +
#      geom_line()+theme_bw()+ggtitle("changepoint") + theme(plot.title = element_text(hjust = 0.5))
#    list<-list(Net_Degree_test,MA2test,MStd2test,MA2MStd2test,BDIquant,t,p)
#    names(list)<-c("BDM","movmean","movsd","BDI","BDIquant","test","figure")
#    return(list)
#  }

## ----fig.width=6, fig.height=4------------------------------------------------
set.seed(1234)
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
  a1<-changepoint(x,8,1)
  a1
  a2<-changepoint(x,6,2)
  a2

## ----fig.width=6, fig.height=4------------------------------------------------
set.seed(1)
y<-rep(0,5000)
epi<-rnorm(5000,0,0.5)
deltay<-c(rep(0,500),rep(1,500),rep(3,500),rep(6,500),rep(10,500),rep(15,500),rep(21,500),rep(28,500),rep(36,500),rep(45,500))
for(i in 3:5000){
  y[i]<-0.6*y[i-1]-0.5*y[i-2]+epi[i]+deltay[i]
}
a<-changepoint(y,15,1)
a$figure

