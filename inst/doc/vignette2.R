## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----cars---------------------------------------------------------------------
plot(seq(0,2*pi,length=200),sin(seq(0,2*pi,length=200)),type='l',col='red',xlab='x',ylab='y')
lines(seq(0,2*pi,length=200),cos(seq(0,2*pi,length=200)),col='blue')
text(2.4,0.7,"sin(x)")
text(5.5,0.7,"cos(x)")
legend(0,-0.5,col=c("red","blue"),lty=c(1,1),lwd=c(1,1),legend=c("sin","cos"))

## -----------------------------------------------------------------------------
A<-matrix(1:9,nrow=3,ncol=3)
B<-diag(1:3)
rbind(solve(B),A)

## -----------------------------------------------------------------------------
C<-crossprod(A,B)
C
C[which(C==min(C,na.rm=T))]
which(C==min(C,na.rm=T),arr.ind=T)

## -----------------------------------------------------------------------------
library(ggplot2)

set.seed(1)
data <- data.frame(
  x = rnorm(100),
  y = rnorm(100)
)

plot <- ggplot(data, aes(x = x, y = y)) +
  geom_point() +                      
  labs(title = "简单散点图",            
       x = "X 轴",                    
       y = "Y 轴")        

print(plot)

## -----------------------------------------------------------------------------
rayleigh<-function(n,sigma){
  u<-runif(n)
  x<-sqrt(-2*sigma^2*log(1-u))
  return(x)
}
ray<-rayleigh(10000,1)
hist(ray,prob="T",breaks=20,main="Histogram of Rayleigh(1)")
y <- seq(0, 4, 0.01)
lines(y, y*exp(-y^2/2))
text(2.5,0.5,expression(f(x)==x*e^{-x^2/2}))
hist(rayleigh(10000,2),prob="T",breaks=20,main="Histogram of Rayleigh(2)")
y <- seq(0, 8, 0.01)
lines(y, y*exp(-y^2/8)/4)
text(5,0.25,expression(f(x)=={(x/4)}*e^{-x^2/8}))

## -----------------------------------------------------------------------------
normal_mixture<-function(n,p){
  x1<-rnorm(n,0,1)
  x2<-rnorm(n,3,1)
  z<-rep(0,n)
  u<-runif(n)
  for(i in 1:n){
    if(u[i]<p){
      z[i]<-x1[i]
    }else{
      z[i]<-x2[i]
    }
  }
  return(z)
}
x<-normal_mixture(1000,0.75)
hist(x,breaks = 20,main="Histogram of normal mixture of p1=0.75")
hist(normal_mixture(1000,0.5),breaks = 30,main="Histogram of normal mixture of p1=0.5")


## -----------------------------------------------------------------------------
pg<-function(t,lambda,a,b,n=1000){
  x<-rep(0,n)
  for(i in 1:n){
    nt<-sum(rpois(t,lambda))
    x[i]<-sum(rgamma(nt,a,b))
  }
  return(x)
}
z<-pg(t=10,lambda=1,a=1,b=2,10000)
hist(z,main="compound Poisson(1)–Gamma process ")
mean(z)
var(z)

## -----------------------------------------------------------------------------
set.seed(12345)
beta33<-function(num){
  m <- 1e4
  x <- runif(m, min=0, num)
  theta.hat <- mean(30*x^2*(1-x)^2)*num
  return(theta.hat)
}
x<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
data<-matrix(rep(0,18),2)
for(i in 1:length(x)){
  data[1,i]<-beta33(x[i])
  data[2,i]<-pbeta(x[i],3,3)
}
round(data,4)

## -----------------------------------------------------------------------------
rayleigh1<- function(n, sigma) {
  u<-runif(n)
  x<-sigma*sqrt(-2*log(1-u))
  x1<-sigma*sqrt(-2*log(u))
  return((x+x1)/2)
}
rayleigh2<-function(n,sigma) {
  X1<-sigma*sqrt(-2*log(runif(n)))
  X2<-sigma*sqrt(-2*log(runif(n)))
  
  return((X1+X2)/2)
}
n <- 10000 
sigma<-1 
#对偶变量
s1<-rayleigh1(n, sigma)
#独立变量
s2<-rayleigh2(n, sigma)
var1<- var(s1)
var1
var2<- var(s2)
var2
(var2-var1)/var2


## -----------------------------------------------------------------------------
set.seed(12345)
g<-function(x){
  (x^2/sqrt(2*pi))*exp(-x^2/2)
}
f1_sample<-function(n){
  u <- runif(n)  
  samples <- qnorm(u * (1 - pnorm(1)) + pnorm(1))#逆变换
  g_values<-g(samples)
  f1_values<-dnorm(samples)/(1-pnorm(1))
  mean<-mean(g_values/f1_values)
  variance<-var(g_values/f1_values)
  return(list(mean,variance))
}
f2_sample<-function(n){
  u <- runif(n)
  samples <- qchisq(u * (1 - pchisq(1, df = 4)) + pchisq(1, df = 4), df = 4)#逆变换
  g_values <- g(samples)
  f2_values <-dchisq(samples,df = 4)/(1-pchisq(1,df = 4))
  mean<-mean(g_values/f2_values)
  variance <- var(g_values/f2_values)
  return(list(mean,variance))
}
n <- 10000
var_f1<-f1_sample(n)
var_f2<-f2_sample(n)
var_f1
var_f2

## -----------------------------------------------------------------------------
library(ggplot2)
n<- c(10^4, 2*10^4, 4*10^4, 6*10^4, 8*10^4)
tn<- n* log(n)
an<- numeric(5)
set.seed(12345) 
#快速排序
quick_sort <- function(n) {
  if (length(n) < 2) {
    return(n)
  } else {
    x <- n[1]
    left<- n[n<x]
    right<- n[n >= x][-1] 
    return(c(quick_sort(left), x, quick_sort(right)))
  }
}

for (i in 1:length(n)) {
  n0<- n[i]
  #100次模拟
  times<- numeric(100)  
  for (j in 1:100) {
    n1 <- sample(1:n0, n0)  
    start<- Sys.time()  
    quick_sort(n1)  
    end<- Sys.time() 
    times[j] <- as.numeric(difftime(end, start, units = "secs"))  
  }
  an[i] <- mean(times)  
}

data <- data.frame(tn = tn, an = an)
model <- lm(an ~ tn, data = data)
summary(model)
ggplot(data, aes(x = tn, y = an)) +
  geom_point(color = 'blue') +  
  geom_smooth(method = "lm", col = "red") + 
  labs(title = "Regression of an on tn",
       x = "tn",
       y = "an")


## -----------------------------------------------------------------------------
set.seed(12345)
n<-1000
skewness<-function(x) {
  m3<-mean((x-mean(x))^3)
  m2<-mean((x-mean(x))^2)
  return(m3/(m2^(3/2)))
}
b1_values<-replicate(n,{
  sn<-rnorm(n)  
  b1<-skewness(sn)
})
qmc<-quantile(b1_values,probs=c(0.025, 0.05, 0.95, 0.975))
qthe<-qnorm(c(0.025, 0.05, 0.95, 0.975),mean=0,sd=sqrt(6/n))
q<-data.frame(
  quantile=c(0.025, 0.05, 0.95, 0.975),
  monteCarlo=qmc,
  theoretical=qthe
)
print(q)


## -----------------------------------------------------------------------------
q_sd<-function(q){
  return(sqrt(q*(1-q)/n/dnorm(qnorm(q,0,sqrt(6/n)),0,sqrt(6/n))))
}
se_mc<-q_sd(c(0.025, 0.05, 0.95, 0.975))
print(se_mc)

## -----------------------------------------------------------------------------
set.seed(12345)
library(MASS)
n <- 100
mu <- c(0, 0) 
Sigma <- matrix(c(1, 0.2, 0.2, 1), 2, 2) 
data_normal <- mvrnorm(n, mu = mu, Sigma = Sigma)
x <- data_normal[, 1]
y <- data_normal[, 2]
cor_pearson <- cor.test(x, y, method = "pearson")
cor_spearman <- cor.test(x, y, method = "spearman")
cor_kendall <- cor.test(x, y, method = "kendall")
cat("双变量正态分布下的相关性检验:\n")
print(cor_pearson)
print(cor_spearman)
print(cor_kendall)
calculate_power <- function(n, rho, alpha = 0.05, method = c("pearson", "spearman", "kendall"), reps = 1000) {
  reject_count <- 0
  for (i in 1:reps) {
    Sigma <- matrix(c(1, rho, rho, 1), 2, 2) 
    data <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
    x <- data[, 1]
    y <- data[, 2]
    test_result <- cor.test(x, y, method = method)
    if (test_result$p.value < alpha) {
      reject_count <- reject_count + 1
    }
  }
  return(reject_count / reps)
} 
rho <- 0.3
alpha <- 0.05  
reps <- 1000 
power_pearson <- calculate_power(n, rho, alpha, method = "pearson", reps = reps)
power_spearman <- calculate_power(n, rho, alpha, method = "spearman", reps = reps)
power_kendall <- calculate_power(n, rho, alpha, method = "kendall", reps = reps)
cat("皮尔逊检验的功效:", power_pearson, "\n")
cat("斯皮尔曼检验的功效:", power_spearman, "\n")
cat("肯德尔检验的功效:", power_kendall, "\n")

## -----------------------------------------------------------------------------
#反例
n<-100
x_alt <- rnorm(n)
y_alt <- x_alt^2 + rnorm(n)
cor_pearson_alt <- cor.test(x_alt, y_alt, method = "pearson")
cor_spearman_alt <- cor.test(x_alt, y_alt, method = "spearman")
cor_kendall_alt <- cor.test(x_alt, y_alt, method = "kendall")
cat("\n双变量t分布的相关性检验:\n")
print(cor_pearson_alt)
print(cor_spearman_alt)
print(cor_kendall_alt)
calculate_power1 <- function(n,rho,df,alpha = 0.05, method = c("pearson", "spearman", "kendall"), reps = 1000) {
  reject_count <- 0
  for (i in 1:reps) {
    Sigma<-matrix(c(1,rho,rho,1),2,2)
    data_t <- mvtnorm::rmvt(n, sigma = Sigma, df = df)
    x_alt <- data_t[, 1]
    y_alt <- data_t[, 2]
    test_result <- cor.test(x_alt, y_alt, method = method)
    if (test_result$p.value < alpha) {
      reject_count <- reject_count + 1
    }
  }
  return(reject_count / reps)
} 
alpha <- 0.05  
reps <- 1000 
rho<-0.3
df<-2
power_pearson <- calculate_power1(n,rho,df,alpha,method = "pearson", reps = reps)
power_spearman <- calculate_power1(n,rho,df,alpha,method = "spearman", reps = reps)
power_kendall <- calculate_power1(n,rho,df,alpha,method = "kendall", reps = reps)
cat("皮尔逊检验的功效:", power_pearson, "\n")
cat("斯皮尔曼检验的功效:", power_spearman, "\n")
cat("肯德尔检验的功效:", power_kendall, "\n")

## -----------------------------------------------------------------------------
d1<- c(-2.961, 0.478, -0.391, -0.869, -0.460, -0.937, 0.779, -1.409, 0.027, -1.569)
d2<- c(1.608, 1.009,  0.878,  1.600, -0.263,  0.680, 2.280,  2.390, 1.793,  8.091,1.468)
B<-1e4
set.seed(12345)
thetastar<-numeric(B)
theta1<-mean(d1)
theta2<-mean(d2)
mean<-(theta1-theta2)
for(b in 1:B){
 xstar1<-sample(d1,length(d1),replace=TRUE)
 xstar2<-sample(d2,length(d2),replace=TRUE)
 thetastar[b]<-mean(xstar1)-mean(xstar2)
 }
round(c(theta.boot=mean(thetastar),theta=mean,bias=mean(thetastar)-mean,se.boot=sd(thetastar),se.samp=sqrt(var(d1)/(length(d1))+var(d2)/(length(d2)))),5)

## -----------------------------------------------------------------------------
set.seed(12345)  
N<-1000      
null_hypotheses<-950
alt_hypotheses<-50
m<-10000
alpha<-0.1
calculate_metrics<-function(p_values){
  reject_null<-p_values<alpha
  FWER<-mean(reject_null[1:null_hypotheses]) 
  TPR<-mean(reject_null[(null_hypotheses+1):N]) 
  FDR<-mean(reject_null)  
  return(c(FWER,FDR,TPR))
}
results<-matrix(0,nrow=m,ncol=6)
for(i in 1:m){
  p_values<-c(runif(null_hypotheses),rbeta(alt_hypotheses,0.1,1))
  bonferroni_adjusted<-p.adjust(p_values,method="bonferroni")
  bonferroni_metrics<-calculate_metrics(bonferroni_adjusted)
  bh_adjusted<-p.adjust(p_values,method="BH")
  bh_metrics<-calculate_metrics(bh_adjusted)
  results[i, ]<-c(bonferroni_metrics,bh_metrics)
}
colMeans(results)

## -----------------------------------------------------------------------------
set.seed(12345) 
library(boot)
times<-c(3,5,7,18,43,85,91,98,100,130,230,487)
mle<-1/mean(times)
bootstrap_lambda<-function(data,i) {
  sample<-data[i]
  return(1/mean(sample))
}
results<-boot(data=times,statistic=bootstrap_lambda,R=1000)
bias<-mean(results$t)-mle
std<-sd(results$t)
mle
bias
std

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
alpha<-0.05
bootstrap_lambda<-function(data,i) {
  sample<-data[i]
  return(mean(sample))
}
results<-boot(data=times,statistic=bootstrap_lambda,R=1000)
mean_boot<-mean(results$t)
std<-sd(results$t)
z<-qnorm(1-alpha/2)
normal<-c(mean_boot-z*std,mean_boot+z*std)
basic<-c(2*mean_boot-quantile(results$t,1-alpha/2),2*mean_boot-quantile(results$t,alpha/2))
percentile<-quantile(results$t, c(alpha/2,1-alpha/2))
bca<-boot.ci(results, type = "bca")$bca[4:5]
normal
basic
percentile
bca

## -----------------------------------------------------------------------------
library(bootstrap)
sigma<-cov(scor)
lambda<-eigen(sigma)$val
lambda1<-lambda[which.max(lambda)]
theta<-lambda1/sum(lambda)
n<-nrow(scor)
theta.jack<-numeric(n)
for(i in 1:n){
  sigma<-cov(scor[-i,])
  lambda<-eigen(sigma)$val
  theta.jack[i]<-lambda[which.max(lambda)]/sum(lambda)
}
bias<-(n-1)*(mean(theta.jack)-theta)
se<-sqrt((n-1)*mean(theta.jack-mean(theta.jack)^2))
bias
se

## -----------------------------------------------------------------------------
 library(DAAG)
 n <- length(ironslag$magnetic) 
 e1 <- e2 <- e3 <- e4<-e5 <- numeric(n)
 for (k in 1:n) {
 y <- ironslag$magnetic[-k]
 x <- ironslag$chemical[-k]
 J1 <- lm(y ~ x)
 yhat1 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[k]
 e1[k] <- ironslag$magnetic[k]- yhat1
 J2 <- lm(y ~ x + I(x^2))
 yhat2 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[k] +
 J2$coef[3] * ironslag$chemical[k]^2
 e2[k] <- ironslag$magnetic[k]- yhat2
 J3 <- lm(log(y) ~ x)
 logyhat3 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[k]
 yhat3 <- exp(logyhat3)
 e3[k] <- ironslag$magnetic[k]- yhat3
 J4 <- lm(log(y) ~ log(x))
 logyhat4 <- J4$coef[1] + J4$coef[2] * log(ironslag$chemical[k])
 yhat4 <- exp(logyhat4)
 e4[k] <- ironslag$magnetic[k]- yhat4
 J5<-lm(y~x+I(x^2)+I(x^3))
 yhat5<-J5$coef[1]+J5$coef[2]*ironslag$chemical[k]+J5$coef[3]*ironslag$chemical[k]^2+J5$coef[4]*ironslag$chemical[k]^3
 e5[k] <- ironslag$magnetic[k]- yhat5
 }
 c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2),mean(e5^2))
 lm(ironslag$magnetic ~ ironslag$chemical + I(ironslag$chemical^2))

## -----------------------------------------------------------------------------
library(DAAG)
 n <- length(ironslag$magnetic) 
 r1 <- r2 <- r3 <- r4<-r5 <- numeric(n)
 for (k in 1:n) {
 y <- ironslag$magnetic[-k]
 x <- ironslag$chemical[-k]
 J1 <- lm(y ~ x)
 r1[k] <- summary(J1)$adj.r.squared
 J2 <- lm(y ~ x + I(x^2))
 r2[k] <- summary(J2)$adj.r.squared
 J3 <- lm(log(y) ~ x)
 r3[k] <- summary(J3)$adj.r.squared
 J4 <- lm(log(y) ~ log(x))
 r4[k] <- summary(J4)$adj.r.squared
 J5<-lm(y~x+I(x^2)+I(x^3))
 r5[k] <- summary(J5)$adj.r.squared
 }
 c(mean(r1^2), mean(r2^2), mean(r3^2), mean(r4^2),mean(r5^2))

## -----------------------------------------------------------------------------
library(twosamples)
 x <- as.vector(chickwts$weight[chickwts$feed == "soybean"])
 y <- as.vector(chickwts$weight[chickwts$feed == "linseed"])
 R <- 999
 z <- c(x, y)
 K <- 1:26
 reps <- numeric(R)
 t0 <- cvm_test(x, y)[1]
 for (i in 1:R) {
 k <- sample(K, size = 14, replace = FALSE)
 x1 <- z[k]
 y1 <- z[-k]
 reps[i] <- cvm_test(x1, y1)[1]
 }
 p <- mean(c(t0, reps) >= t0)
 p

## -----------------------------------------------------------------------------
x <- as.vector(chickwts$weight[chickwts$feed == "sunflower"])
 y <- as.vector(chickwts$weight[chickwts$feed == "linseed"])
 R <- 999
 z <- c(x, y)
 K <- 1:24
 reps <- numeric(R)
 t0 <- cor(x, y,method="spearman")
 for (i in 1:R) {
 k <- sample(K, size = 12, replace = FALSE)
 x1 <- z[k]
 y1 <- z[-k]
 reps[i] <- cor(x1, y1,method="spearman")
 }
 p <- mean(c(t0, reps) >= t0)
 p
 cor.test(x, y,method="spearman")$p.value

## -----------------------------------------------------------------------------
set.seed(1234)
library(ggplot2)
metropolis_hastings <- function(n) {
  x<-numeric(n)
  x[1]<-rnorm(1,0,sd)  
  for (i in 2:n) {
    y <- rnorm(1,mean=x[i-1],sd)
    ratio<-dcauchy(y)*dnorm(x[i-1],mean=y,sd)/(dcauchy(x[i-1])*dnorm(y,mean=x[i-1],sd))
    if (runif(1)>ratio) {
      x[i]<-x[i-1]
    } else {
      x[i]<-y
    }
  }
  return(x)
}
gelman_rubin <- function(chains) {
  m <- ncol(chains)  
  n <- nrow(chains)
  means <- colMeans(chains)
  overall_mean <- mean(means)
  B <- n * var(means)  #序列间方差
  W <- mean(apply(chains, 2, var))  #序列内方差
  var_hat <- (n - 1) / n * W + B / n
  R_hat <- (var_hat / W)
  return(R_hat)
}
n_samples <- 10000
n_chains <- 4
sd<-1
k<-1
chains <- matrix(NA, n_samples, n_chains)
for (i in 1:n_chains) {
  chains[, i] <- metropolis_hastings(n_samples)
}
R_hat<-numeric(n_samples)
repeat {
  k<-k+1
  if(k>n_samples){
    break
  }
  R_hat[k] <- gelman_rubin(chains[1:k,])
}
plot(R_hat,type='l',ylim=c(0,2))
abline(h=1.2)


## -----------------------------------------------------------------------------
chains <- chains[-(1:1000), ]
generated_deciles <- apply(chains, 2, quantile, probs = seq(0.1, 0.9, by = 0.1))
theoretical_deciles <- qcauchy(seq(0.1, 0.9, by = 0.1))
comparison <- data.frame(
  Decile = seq(0.1, 0.9, by = 0.1),
  Generated = rowMeans(generated_deciles),
  Theoretical = theoretical_deciles
)
print(comparison)
ggplot(comparison, aes(x = Decile)) +
  geom_line(aes(y = Generated, color = "Generated"), linewidth = 1) +
  geom_line(aes(y = Theoretical, color = "Theoretical"), linewidth = 1) +
  labs(title = "Decile Comparison: Generated vs Theoretical Cauchy Distribution",
       x = "Decile", y = "Value") +
  scale_color_manual(values = c("Generated" = "blue", "Theoretical" = "red")) +
  theme_minimal()


## -----------------------------------------------------------------------------
library(ggplot2)
gibbs<- function(n_samples) {
  n <- 10 
  x <- numeric(n_samples)
  y<- numeric(n_samples)
  y[1]<-0.5
  x[1]<-0.5
  for (i in 2:n_samples) {
    x1<-x[i-1]
    y1<-y[i-1]
    x1<- rbinom(1, n, y1)
    x[i]<-x1
    y[i]<- rbeta(1, x1+a,n-x1+ b)
  }
  return(data.frame(x,y))
}
    
gelman_rubin <- function(chains) {
  m <- ncol(chains)  
  n <- nrow(chains)
  means <- colMeans(chains)
  overall_mean <- mean(means)
  B <- n * var(means)  #序列间方差
  W <- mean(apply(chains, 2, var))  #序列内方差
  var_hat <- (n - 1) / n * W + B / n
  R_hat <- (var_hat / W)
  return(R_hat)
}
a <- 2     
b <- 3 
n_samples <- 10000
n_chains <- 4
k<-1
n<-10
chainsx<-chainsy<-matrix(NA, n_samples, n_chains)
for (i in 1:n_chains) {
  chainsx[,i]<-gibbs(n_samples)$x
  chainsy[,i]<-gibbs(n_samples)$y
}
R_hatx<-R_haty<-rep(0,n_samples)
repeat {
  k<-k+1
  if(k>n_samples){
    break
  }
  R_hatx[k] <- gelman_rubin(chainsx[1:k,])
  R_haty[k] <- gelman_rubin(chainsy[1:k,])
}
plot(R_hatx[1:n_samples],type='l')
plot(R_haty[1:n_samples],type='l')
par(mfrow=c(1, 2))
hist(chainsx[,1],breaks = seq(-0.5, n + 0.5, by = 1), probability = TRUE,main = "Histogram of x samples",xlab = "x", col = "lightblue", xlim = c(0, n))
hist(chainsy[,1], breaks = 30, probability = TRUE,main = "Histogram of y samples",xlab = "y", col = "lightgreen", xlim = c(0, 1))

## -----------------------------------------------------------------------------
fk<-function(k,a,d){
  x<-(-1)^k/2^k/factorial(k)
  b<-norm(a,"2")^(2*k+2)/(2*k+1)/(2*k+2)
  c<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  return(x*b*c)
}

## -----------------------------------------------------------------------------
f<-function(a,d){
  sumk<-0
  for(k in 0:50){
    x<-(-1)^k/2^k/factorial(k)
    b<-norm(a,"2")^(2*k+2)/(2*k+1)/(2*k+2)
    c<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
    sumk<-sumk+x*b*c
  }
  return(sumk)
}

## -----------------------------------------------------------------------------
f<-function(a,d){
  sumk<-0
  for(k in 0:100){
    x<-(-1)^k/2^k/factorial(k)
    b<-norm(a,"2")^(2*k+2)/(2*k+1)/(2*k+2)
    c<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
    sumk<-sumk+x*b*c
  }
  return(sumk)
}
a<-c(1,2)
d<-1
f(a,d)

## -----------------------------------------------------------------------------
left_integral<-function(ck1, k){
  integrand<-function(u) (1+u^2/(k-1))^(-k/2)
  integrate(integrand,0,ck1)$value
}
right_integral<-function(ck, k){
  integrand<-function(u) (1 + u^2/k)^(-(k + 1) / 2)
  integrate(integrand,0,ck)$value
}
equation<-function(a, k) {
  ck1<-sqrt(a^2*(k-1)/(k-a^2))
  ck<-sqrt(a^2*k/(k+1-a^2))
  left<-(2*exp(lgamma((k)/2)-lgamma((k-1)/2))/sqrt(pi*(k-1)))*left_integral(ck1, k)
  right<-(2*exp(lgamma((k+1)/2)-lgamma(k/2))/sqrt(pi*k))*right_integral(ck, k)
  return(left-right)
}
solve_for_a <- function(k) {
  solution <- uniroot(equation,interval=c(0.1,sqrt(k)-1),k=k)
  return(solution$root)
}
k<-10
a_solution<-solve_for_a(k)
a_solution

## -----------------------------------------------------------------------------
Y<-c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau<-1
n<-length(Y)
epsilon<-0.0001
lambda<-1
repeat {
  lambda_old<-lambda
  Q_t<-sum(Y[Y<tau])+sum((tau+lambda)*(Y==tau))
  lambda<-Q_t/n
  if (abs(lambda-lambda_old)< epsilon){
    break
  }
}
lambda


## -----------------------------------------------------------------------------
lambdamle<-(sum(Y))/sum(Y<1)
lambdamle

## -----------------------------------------------------------------------------
library(lpSolve)
ob<-c(4, 2, 9)
con<-matrix(c(
  2, 1, 1,  
  1, -1, 3  
),nrow=2, byrow=TRUE)
rhs<-c(2, 3)
dir<-c("<=","<=")
sol<-lp("min",ob,con,dir,rhs)
sol$solution 
sol$objval   

## -----------------------------------------------------------------------------
formulas<-list(
  mpg~disp,
  mpg~I(1/disp),
  mpg~disp+wt,
  mpg~I(1/disp)+wt
)
forloop3<-list()
for (i in seq_along(formulas)) {
  forloop3[[i]]<-lm(formulas[[i]],data=mtcars)
}
forloop3
lapply3<- lapply(formulas, function(f) lm(f, data = mtcars))
lapply3

## -----------------------------------------------------------------------------
set.seed(12345) 
bootstraps<-lapply(1:10,function(i){
  rows<-sample(1:nrow(mtcars),rep=TRUE)
  mtcars[rows, ]
})
forloop4<-list()
for(i in seq_along(bootstraps)){
  forloop4[[i]]<-lm(mpg~disp,data=bootstraps[[i]])
}
forloop4
fit<-function(data){
  lm(mpg~disp,data=data)
}
lapply4<-lapply(bootstraps,fit)
lapply4

## -----------------------------------------------------------------------------
rsq<-function(mod) summary(mod)$r.squared
rsqfor3<-numeric(length(forloop3))
for(i in seq_along(forloop3)){
  rsqfor3[i]<-rsq(forloop3[[i]])
}
rsqfor3
rsqlapply3<-sapply(lapply3, rsq)
rsqlapply3

## -----------------------------------------------------------------------------
rsq<-function(mod) summary(mod)$r.squared
rsqfor4<-numeric(length(forloop4))
for(i in seq_along(forloop4)){
  rsqfor4[i]<-rsq(forloop4[[i]])
}
rsqfor4
rsqlapply4<-sapply(lapply4, rsq)
rsqlapply4

## -----------------------------------------------------------------------------
trials<-replicate(
  100,
  t.test(rpois(10,10),rpois(7,10)),
  simplify=FALSE
)
p_values<-sapply(trials, function(t) t$p.value)
p_values
p_values<-sapply(trials, `[[`, "p.value")
p_values

## -----------------------------------------------------------------------------
lapply_variant<-function(FUN,...,FUN.VALUE){
  inputs<-list(...) 
  Map(FUN, ...)|>vapply(identity,FUN.VALUE=FUN.VALUE)
}
x<-1:10
y<-11:20
result<-lapply_variant(function(a,b) a+b,x,y,FUN.VALUE=numeric(1))
print(result) 
result_matrix<-lapply_variant(function(a,b) c(sum=a+b,product=a*b),x,y,FUN.VALUE=numeric(2))
print(result_matrix)

## ----warning=FALSE------------------------------------------------------------
fast_chisq1<-function(x, y){
  con_table<-table(x,y)
  row<-rowSums(con_table)
  col<-colSums(con_table)
  grand<-sum(con_table)
  expected<-outer(row,col)/grand
  observed<-as.numeric(con_table)
  expected<-as.numeric(expected)
  chisq_stat<-sum((observed-expected)^2/expected)
  return(chisq_stat)
}
fast_chisq2<-function(x,y){
  con_table<-table(x,y)
  E<-outer(rowSums(con_table),colSums(con_table))/sum(con_table)
  O<-con_table
  chisq_stat<-sum((O-E)^2/E)
  return(chisq_stat)
}
x<-c(1,2,2,3,1,3,2)
y<-c(2,2,3,1,1,3,2)
chisq_fast1<-fast_chisq1(x,y)
print(chisq_fast1)
chisq_fast2<-fast_chisq2(x,y)
print(chisq_fast2)
chisq_original<-chisq.test(table(x,y))$statistic
print(chisq_original)

## ----warning=FALSE------------------------------------------------------------
fast_table<-function(x, y){
  x_levels<-sort(unique(x))
  y_levels<-sort(unique(y))
  result<-matrix(0L,nrow=length(x_levels),ncol=length(y_levels),dimnames=list(x_levels,y_levels))
  for(i in seq_along(x)){
    result[as.character(x[i]),as.character(y[i])]<-result[as.character(x[i]),as.character(y[i])]+1
  }
  return(result)
}
fast_chisq<-function(x,y){
  con_table<-fast_table(x, y)
  row<-rowSums(con_table)
  col<-colSums(con_table)
  grand<-sum(con_table)
  expected<-outer(row,col)/grand
  observed<-as.numeric(con_table)
  expected<-as.numeric(expected)
  chisq_stat<-sum((observed-expected)^2/expected)
  return(chisq_stat)
}
x<-as.integer(c(1,2,2,3,1,3,2))
y<-as.integer(c(2,2,3,1,1,3,2))
table_fast<-fast_table(x,y)
table_fast
chisq_fast<-fast_chisq(x, y)
chisq_fast
table_original<-table(x, y)
chisq_original<-chisq.test(table(x, y))$statistic
chisq_original

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction("
NumericMatrix gibbs_sampler(int n,double a,double b,int N,int init_x,double init_y) {
  NumericMatrix samples(N,2);
  int x=init_x;
  double y=init_y;
  for(int i=0;i<N;i++){
    x=R::rbinom(n,y);
    y=R::rbeta(x+a,n-x+b);
    samples(i,0)=x;
    samples(i,1)=y;
  }
  return samples;
}
")
n<-10
a<-2
b<-3
N<-10000
init_x<-0.5
init_y<-0.5
samples<-gibbs_sampler(n,a,b,N,init_x,init_y)
par(mfrow=c(1,2))
plot(samples[,2],type="l",xlab="Iteration",ylab="y")
hist(samples[,2],breaks=30,main="y",xlab="y")

## -----------------------------------------------------------------------------
gibbs<-function(N){
  n<-10 
  x<-numeric(N)
  y<-numeric(N)
  y[1]<-0.5
  x[1]<-0.5
  for(i in 2:N){
    x1<-x[i-1]
    y1<-y[i-1]
    x1<-rbinom(1,n,y1)
    x[i]<-x1
    y[i]<-rbeta(1,x1+a,n-x1+b)
  }
  return(data.frame(x,y))
}
n<-10
a<-2
b<-3
N<-10000
init_x<-0.5
init_y<-0.5
samplescpp<-gibbs_sampler(n,a,b,N,init_x,init_y)  
samplesR<-gibbs(N)
qqplot(samplescpp[,1],samplesR[,1])
qqplot(samplescpp[,2],samplesR[,2])

## -----------------------------------------------------------------------------
library(microbenchmark)
ts<-microbenchmark(gibbs_sampler(n,a,b,N,init_x,init_y),gibbs(N))
summary(ts)

