## -----------------------------------------------------------------------------
SCAD_cv <- function(X, y, A) {

  #to define SCAD penalty derivative
  pen_der<-function(theta,a = 3.7,lambda){
    if(abs(theta) > lambda){
      if(a * lambda > theta){
        return((a * lambda - theta)/(a - 1))
      }else{
        return(0)
      }
    }else{
      return(lambda)
    }
  }
  
  #Initialization coefficient
  b_ls <- stats::lm(y~X)$coefficients
  #the maximum iterations  
  max_iteration<-1e4
  #A column 1 is added before X to facilitate the processing of intercept items
  inter_y <- matrix(1,nrow = length(y),ncol = 1)
  X <- cbind(inter_y,X)
  
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  CV_score <- rep(0,length(A))
  
  for(i in 1:length(A)){
  b0 <- b_ls
  ll <- A[i]
  iteration<-0
  while(iteration < max_iteration){
    for(j in 1:length(b0)){
      if(abs(b0[j]) < 1e-06){
        next()
      }else{
        nominator<-sum(- X[,j] * (y - X %*% b0)) + nrow(X)*b0[j]*pen_der(theta = b0[j],lambda =ll)/abs(b0[j])
        denominator<- sum(X[,j]^2) + nrow(X)*pen_der(theta = b0[j],lambda =ll)/abs(b0[j])
        b0[j]<-b0[j] - nominator/denominator
        if(abs(b0[j]) < 1e-06){
          b0[j] <- 0
        }
      }
    }
    iteration <- iteration+1 
  }
  
  y_hat <- X %*% b0
  s = 0
  for(k in 1:length(y)){
    s = s + ((y[k]-y_hat[k])/(1-H[k,k]))^2
  }
  CV_score[i] <- s/length(y)
  }
  
  for(i in 1:length(A)){
    if(CV_score[length(A)-i+1] < 1){
      index=(length(A)-i+1)
      break
    }
  }
  return( A[index])
}

## -----------------------------------------------------------------------------
SCAD_Function <- function(X, y) {
  
  #to find the best lambda by cv
  b0<- stats::lm(y~X)$coefficients
  qq <- max(abs(b0))
  if(qq > 50){
    A <- seq(0,qq,1)
  }else{
    A <- seq(0,qq,0.1)
  }
  cv_lambda <- SCAD_cv(X,y,A)
  
  #to define SCAD penalty
  pen_fun<-function(theta,lambda = cv_lambda){
    p_lambda<-sapply(theta, function(x){
      if(abs(x)< lambda){
        return(lambda^2 - (abs(x) - lambda)^2)
      }else{
        return(lambda^2)
      }
    }
    )
    return(p_lambda)
  }
  
  #to define SCAD penalty derivative
  pen_der<-function(theta,a = 3.7,lambda = cv_lambda){
    if(abs(theta) > lambda){
      if(a * lambda > theta){
        return((a * lambda - theta)/(a - 1))
      }else{
        return(0)
      }
    }else{
      return(lambda)
    }
  }
  
  #Initialization coefficient
  b0<- stats::lm(y~X)$coefficients
  #Iteration counter
  iteration<-0
  #the maximum iterations  
  max_iteration<-100000
  #A column 1 is added before X to facilitate the processing of intercept items
  inter_y <- matrix(1,nrow = length(y),ncol = 1)
  X <- cbind(inter_y,X)
  
  #SCAD iterative process
  while(iteration < max_iteration){
    for(j in 1:length(b0)){
      if(abs(b0[j]) < 1e-06){
        next()
      }else{
        #nominator
        nominator<-sum(- X[,j] * (y - X %*% b0)) + nrow(X)*b0[j]*pen_der(theta = b0[j])/abs(b0[j])
        #denominator__the second derivative
        denominator<- sum(X[,j]^2) + nrow(X)*pen_der(theta = b0[j])/abs(b0[j])
        #Iterate b0[j]
        b0[j]<-b0[j] - nominator/denominator
        #If b0[j] is small, make it equal to 0.
        if(abs(b0[j]) < 1e-06){
          b0[j] <- 0
        }
      }
    }
    iteration <- iteration+1 
  }
  return(b0)
}

## -----------------------------------------------------------------------------
library(MASS)
library(ncvreg)
library(stats)

Simu_Multi_Norm<-function(x_len, sd = 1, pho = 0.5){
  V <- matrix(data = NA, nrow = x_len, ncol = x_len)
  for(i in 1:x_len){ 
    for(j in 1:x_len){ 
      V[i,j] <- pho^abs(i-j)
    }
  }
  V<-(sd^2) * V
  return(V)
}

set.seed(123)
X<-mvrnorm(n = 60, mu = rep(0,8), Simu_Multi_Norm(x_len = 8,sd  = 1, pho = 0.5))
beta<-c(1,0.5,0.02,0,2,0,10,0.05)
y<- X %*% beta + rnorm(n = 60, 0, 1)

A <- seq(0,2,0.1)
SCAD_cv(X,y,A)

B <- cv.ncvreg(X,y,penalty='SCAD')
B$lambda.min

SCAD_Function(X,y)


## -----------------------------------------------------------------------------
A<-c(172,170,180,175,188,168,177,183,173)
B<-c(65,80,75,60,85,70,68,77,63)
plot(A,B,xlab="Height/cm",ylab="Weight/kg",xlim=c(160,200),ylim=c(55,90),main="Height and weight of 9 boys in a class")

## -----------------------------------------------------------------------------
a<-matrix(runif(15),5,3)
dimnames(a)<-list(NULL,c("Group one","Group two","Group three"))
knitr::kable(head(a))

## -----------------------------------------------------------------------------

n <- 1000; qq <- 1
X <- sqrt(-2*(qq^2)*log(1-runif(n)))
hist(X, prob = TRUE, main = "f(x),when qq=1")
Y <- seq( 0, 5, .01)
lines(Y, (Y/(qq^2))*exp(-(Y^2)/(2*qq^2)))

n <- 1000; qq <- 2
X <- sqrt(-2*(qq^2)*log(1-runif(n)))
hist(X, prob = TRUE, main = "f(x),when qq=2")
Y <- seq( 0, 10, .01)
lines(Y, (Y/(qq^2))*exp(-(Y^2)/(2*qq^2)))

n <- 1000; qq=3
X <- sqrt(-2*(qq^2)*log(1-runif(n)))
hist(X, prob = TRUE, main = "f(x),when qq=3")
Y <- seq( 0, 15, .01)
lines(Y, (Y/(qq^2))*exp(-(Y^2)/(2*qq^2)))

n <- 1000; qq=4
X <- sqrt(-2*(qq^2)*log(1-runif(n)))
hist(X, prob = TRUE, main = "f(x),when qq=4")
Y <- seq( 0, 20, .01)
lines(Y, (Y/(qq^2))*exp(-(Y^2)/(2*qq^2)))

## ----echo=FALSE---------------------------------------------------------------
n <- 1000 ; p1 <- 0.75
Z1 <- rnorm(n,0,1)
Z2 <- rnorm(n,3,1)
A <- runif(n)
k <- as.integer(A < p1)
X <- k * Z1 + (1-k) * Z2
hist(X, prob = TRUE,main=expression(p[1]==0.75))

## -----------------------------------------------------------------------------
n <- 1000 ; 
Z1 <- rnorm(n,0,1)
Z2 <- rnorm(n,3,1)

A1 <- runif(n)
k1 <- as.integer(A1 < 0.5)
X1 <- k1 * Z1 + (1-k1) * Z2

A2 <- runif(n)
k2 <- as.integer(A2 < 0.25)
X2 <- k2 * Z1 + (1-k2) * Z2

A3 <- runif(n)
k3 <- as.integer(A3 < 0.9)
X3 <- k3 * Z1 + (1-k3) * Z2

A4 <- runif(n)
k4 <- as.integer(A4 < 0.1)
X4 <- k4 * Z1 + (1-k4) * Z2


hist(X1, prob = TRUE, main=expression(p[1]==0.5))
hist(X2, prob = TRUE, main=expression(p[1]==0.25))
hist(X3, prob = TRUE, main=expression(p[1]==0.9))
hist(X4, prob = TRUE, main=expression(p[1]==0.1))

## ----echo=FALSE---------------------------------------------------------------
n=1000;lambda=2;r=3;beta=4
x <- numeric(n)
for(i in 1:n)
{
  a <- rpois(1,10*lambda)
  y <- rgamma(a,r,beta)
  x[i]=sum(y)
}
mean(x)
var(x)


## ----echo=FALSE---------------------------------------------------------------
n=1000;lambda=5;r=6;beta=10
x <- numeric(n)
for(i in 1:n)
{
  a <- rpois(1,10*lambda)
  y <- rgamma(a,r,beta)
  x[i]=sum(y)
}
mean(x)
var(x)

## -----------------------------------------------------------------------------
n <- 1000;

f <- function(x)
{
  a <- runif(n, min=0, max=x);
  b <- mean(a^2*(1-a)^2)*x
  return(b)
}

c <- numeric(9)
d <- numeric(9)

for(i in 1:9)
{
  c[i]=f(i/10)*30
  d[i]=pbeta(i/10,3,3)
}

e <-c((1:9)/10)
p <- rbind(e,c,d)
rownames(p) <- c("x","estimates","pbeta")
knitr::kable(head(p))


## -----------------------------------------------------------------------------
rf <- function(sigma)
{
  n <- 1000;
  a <- runif(n/2)
  b <- runif(n/2)
  X <- sqrt(-2*(sigma^2)*log(1-a))
  X_anti <- sqrt(-2*(sigma^2)*log(a))

  X1 <- sqrt(-2*(sigma^2)*log(1-a))
  X2 <- sqrt(-2*(sigma^2)*log(1-b))

  e <- var((X+X_anti)/2)
  f <- var((X1+X2)/2)
  
  return((f-e)/f)
}

g <- c(1:10)
h <- numeric(10)
for(i in 1:10)
  h[i] <- rf(i)

k <- rbind(g,h)
k <- round(k,3)
rownames(k) <- c("σ","reduction in var")
knitr::kable(head(k))

## -----------------------------------------------------------------------------
n <- 1000
theta.hat <- se <- numeric(2)
g <- function(x)
{
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x>1)
}

x <- rexp(n,1)    #using f1
fg <- g(x)/exp(-x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

x <- rnorm(n,1.65,1)   #using f2
fg <- g(x)/(exp(-(x-1.65)^2/2)/sqrt(2*pi))
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

c <- rbind(theta.hat,se)
colnames(c) <- paste0('f',1:2)
knitr::kable(c)


## -----------------------------------------------------------------------------
n <- 1000
g <- function(x)
{
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x>1)
}

x <- rnorm(n,1.65,1)   #using f2
fg <- g(x)/(exp(-(x-1.65)^2/2)/sqrt(2*pi))
mean(fg)

## -----------------------------------------------------------------------------
set.seed(1000)
n <- 20;alpha <- 0.05
a <- numeric(1000) ;b <- numeric(1000)
for(i in 1:1000)
{
  x <- rchisq(n,df=2)
  a[i] <- mean(x)-qt(alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
  b[i] <- mean(x)+qt(alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)
}

mean( a>2 & b<2 )

## -----------------------------------------------------------------------------
set.seed(10000)
n <- 20;alpha <- 0.05;u=1

a <- numeric(100) ;b <- numeric(100) ;c <- numeric(100)

for(i in 1:100)
{
  x <- rchisq(n,df=1)
  ttest1 <- t.test(x, mu=u)
  a[i] <- ttest1$p.value
  
  y <- runif(n,0,2)
  ttest2 <- t.test(y,alternative = "two.sided", mu=u)
  b[i] <- ttest2$p.value
  
  z <- rexp(n,1)
  ttest3 <- t.test(z,mu=u)
  c[i] <- ttest3$p.value
}

d <- e <- numeric(3)
d[1] <- mean(a<alpha)
d[2] <- mean(b<alpha)
d[3] <- mean(c<alpha)

e[1] <- sqrt(d[1]*(1-d[1])/100)
e[2] <- sqrt(d[2]*(1-d[2])/100)
e[3] <- sqrt(d[3]*(1-d[3])/100)

p<- rbind(d,e)
rownames(p) <- c("alpha.hat","se.hat")
colnames(p) <- c("卡方(1)","Uniform(0,2)","Exponential(1)")
knitr::kable(p)


## -----------------------------------------------------------------------------
library(MASS)
set.seed(29999)
d <- 2
u <- c(0,0)        
v1 <- c(1,0)
v2 <- c(0,2)
v <- rbind(v1,v2)

n <- c(10,20,30,50,100,500)
cr <- qchisq(0.95 , d*(d+1)*(d+2)/6)



T1 <- function(x){          
  n<-nrow(x)
  x1 <- matrix(0,n,d)      
  for(i in 1:d){
       x1[,i]<-x[,i]-mean(x[,i])
  }
  covbar <- matrix(0,d,d)
  for(i in 1:n){
    covbar <- covbar + outer(x1[i,],x1[i,])
  }
  covbar <- covbar/n 
  nicovbar <- solve(covbar)
  b_id <- 0                
  for(i in 1:n){
    for(j in 1:n){
      b_id <- b_id + (x1[i,]%*%nicovbar%*%x1[j,])^3
    }
  }
  b_id <- b_id / (n*n)
  return(n*b_id/6)
}
Tl <-function(x){
    n<-nrow(x)
    x1<-matrix(0,n,d)
    for(i in 1:d){
       x1[,i]<-x[,i]-mean(x[,i])
    }
    covbar<-t(x1)%*%x1/n
    y<-x1%*%solve(covbar)%*%t(x1)
    z<-sum(colSums(y^{3}))/(n*n)  
    return(n*z/6)
}

a <- numeric(6)
b <- 10
for(i in 1:6){
  s<- numeric(b)
  for(j in 1:b){
    x <- mvrnorm(n[i],u,v)
    s[j] <- as.integer(abs(Tl(x)) >= cr)
  }
  a[i] <- mean(s)
}

q <-rbind(n,a)
rownames(q) <- c("n","estimate")
knitr::kable(head(q))


## -----------------------------------------------------------------------------
library(MASS)
set.seed(19999)

u <- c(0,0)
d = 2
afa = 0.1
a = 30
b = 199
ep <- c(seq(0,0.15,0.1),seq(0.15,1,0.5))
L <- length(ep)
po <- numeric(L)
cr <- qchisq(0.9,d*(d+1)*(d+2)/6)  
x <- matrix(0,a,d)

for(j in 1:L){       
  q <- ep[j]
  s <- numeric(b)
  for(i in 1:b){
    σ <- sample(c(1,64), replace = T , size = a , prob = c(1-q , q) )  
       for(g in 1:a){                                              
         r <- matrix(c(σ[g],0,0,σ[g]),2,2)
         x[g,] <- mvrnorm(1, u, r)
       }
    s[i] <- as.integer(abs(Tl(x))>=cr)
  }
  po[j] <- mean(s)
}

plot(ep,po,type = "b" , xlab = bquote(ep) , ylim = c(0,1))    
abline(h = 0.1, lty=3)
stanerror <- sqrt(po*(1-po)/b)
lines(ep,po+stanerror,lty=3)
lines(ep,po-stanerror,lty=3)

## -----------------------------------------------------------------------------
data(scor,package="bootstrap")
m <- nrow(scor)
n <- ncol(scor)
x <- numeric(n)          
for(i in 1:n){                
  x[i] = mean(scor[,i])
}
y <- matrix(0,m,n)
for(i in 1:n){
  y[,i]=scor[,i]-x[i]
}
z <- t(y)%*%y/(m-1)

a <- eigen(z)$values     
t = a[1]/sum(a) 
t

library(boot)
set.seed(10000)
f <- function(x,i){       
  y <- x[i,]
  a <- eigen(cov(y))$values 
  return(a[1]/sum(a)) 
}

bhat = boot(data = scor, statistic = f, R = 2000)  
bis_bo = mean(bhat$t)-t        
se_bo = sd(bhat$t)

bis_bo
se_bo

## -----------------------------------------------------------------------------
x <- numeric(m)
for(i in 1:m){         
  y <- eigen(cov(scor[-i,]))$values
  x[i] <- y[1]/sum(y)
}
bis_ja = (m-1) * (mean(x)-t)
se_ja = (m-1) * sqrt(var(x)/m)

bis_ja
se_ja

## -----------------------------------------------------------------------------
boot.ci(bhat, type = c("perc","bca"))

## -----------------------------------------------------------------------------
library(boot)
set.seed(30000)

f <- function(x,i){   
  y=x[i]
  ybar <- mean(y)
  p <- mean((y-ybar)^3)
  q <- mean((y-ybar)^2)
  return(p/(q^(3/2)))
}

a1 = b1 = c1 = 0      
a2 = b2 = c2 = 0
a3 = b3 = c3 = 0

for(i in 1:10){      
  x <- rnorm(10)
  y <- 0               
  bhat = boot(data = x, statistic = f, R = 500)  
  ci = boot.ci(bhat, type = c("norm","basic","perc"))
  if( y < ci$norm[2]){
    a1 = a1 + 1;
    c1 = c1 + 1;
  }
  if( y > ci$norm[3]){
    b1 = b1 + 1;
    c1 = c1 + 1;
  }
  if( y < ci$basic[4]){
    a2 = a2 + 1;
    c2 = c2 + 1;
  }
  if( y > ci$basic[5]){
    b2 = b2 + 1;
    c2 = c2 + 1;
  }
  if( y < ci$perc[4]){
    a3 = a3 + 1;
    c3 = c3 + 1;
  }
  if( y > ci$perc[5]){
    b3 = b3 + 1;
    c3 = c3 + 1;
  }
}

A <-c(1-c1/10,a1/10,b1/10)
B <-c(1-c2/10,a2/10,b2/10)
C <-c(1-c3/10,a3/10,b3/10)

Z <- rbind(A,B,C)
rownames(Z) <- c("norm","basic","perc")
colnames(Z) <- c("covskewneserage rate","left miss","right miss")
knitr::kable(head(Z))

## -----------------------------------------------------------------------------
library(boot)
set.seed(30000)

f <- function(x,i){       
  y=x[i]
  ybar <- mean(y)
  p <- mean((y-ybar)^3)
  q <- mean((y-ybar)^2)
  return(p/(q^(3/2)))
}

a1 = b1 = c1 = 0      
a2 = b2 = c2 = 0
a3 = b3 = c3 = 0

for(i in 1:10){  
  x <- rchisq(10,5)
  y <-  sqrt(1.6)               
  bhat = boot(data = x, statistic = f, R = 500)  
  ci = boot.ci(bhat, type = c("norm","basic","perc"))
  if( y < ci$norm[2]){
    a1 = a1 + 1;
    c1 = c1 + 1;
  }
  if( y > ci$norm[3]){
    b1 = b1 + 1;
    c1 = c1 + 1;
  }
  if( y < ci$basic[4]){
    a2 = a2 + 1;
    c2 = c2 + 1;
  }
  if( y > ci$basic[5]){
    b2 = b2 + 1;
    c2 = c2 + 1;
  }
  if( y < ci$perc[4]){
    a3 = a3 + 1;
    c3 = c3 + 1;
  }
  if( y > ci$perc[5]){
    b3 = b3 + 1;
    c3 = c3 + 1;
  }
}

A <-c(1-c1/10,a1/10,b1/10)
B <-c(1-c2/10,a2/10,b2/10)
C <-c(1-c3/10,a3/10,b3/10)

Z <- rbind(A,B,C)
rownames(Z) <- c("norm","basic","perc")
colnames(Z) <- c("coverage rate","left miss","right miss")
knitr::kable(head(Z))

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)

xyCor <- function(Z,ix){
  x <- Z[,1]
  y <- Z[ix,2]
  return(cor(x,y,method = "spearman"))
}

X1 <- rnorm(15,0,1)
Y1 <- X1+1
Z1 <- cbind(X1,Y1)
a1 <- boot(data = Z1, statistic =xyCor, R=999, sim="permutation")
b1 <- c(a1$t0,a1$t)
mean(b1 >= a1$t0)
c1 <- cor.test(X1,Y1)
c1$p.value

## -----------------------------------------------------------------------------
X2 <- rnorm(15,0,1)    
Y2 <- rchisq(15,2)
Z2 <- cbind(X2,Y2)
a2 <- boot(data = Z2, statistic =xyCor, R=999, sim="permutation")
b2 <- c(a2$t0,a2$t)
mean(b2 >= a2$t0)
c2 <- cor.test(X2,Y2)
c2$p.value

## -----------------------------------------------------------------------------
set.seed(12345)
library(RANN)
library(boot)
library(energy)
library(Ball)

alpha=0.1;
N <- 100; k<-3; p<-2;B<-99; 
L1 <- L2 <- 15; L <- L1+L2; LL = c(L1,L2)

NN <- function(X, ix, sizes,k) {
  L1 <- sizes[1]; L2 <- sizes[2]; L <- L1 + L2
  if(is.vector(X)) X <- data.frame(X,0);
  X <- X[ix, ];
  Y <- nn2(data=X, k=k+1) 
  part1 <- Y$nn.idx[1:L1,-1]
  part2 <- Y$nn.idx[(L1+1):L,-1]
  i <- sum(part1 <= L1); j <- sum(part2 > L1)
  (i + j) / (k * L)
}

NNtest <- function(X,sizes,k){
  b <- boot(data=X,statistic=NN,R=B,sim = "permutation", sizes = sizes,k=k)
  T <- c(b$t0,b$t)
  pv <- mean(T>=T[1])
  list(statistic=T[1],p.value=pv)
}

p1 <- matrix(NA,N,3)
for(i in 1:N){
  A1 <- matrix(rnorm(L1*p),ncol=p);
  B1 <- cbind(rnorm(L2,0,2),rnorm(L2,0,2));
  C1 <- rbind(A1,B1)
  p1[i,1] <- NNtest(C1,LL,k)$p.value
  p1[i,2] <- eqdist.etest(C1,sizes=LL,R=B)$p.value
  p1[i,3] <- bd.test(x=A1,y=B1,num.permutations = B,seed=i*12345)$p.value
}
powercomparison1 <- colMeans(p1<alpha)

p2 <- matrix(NA,N,3)
for(i in 1:N){
  A2 <- cbind(rnorm(L1,0,2),rnorm(L1,0,2));
  B2 <- cbind(rnorm(L2,1.5,3),rnorm(L2,1.5,3));
  C2 <- rbind(A2,B2)
  p2[i,1] <- NNtest(C2,LL,k)$p.value
  p2[i,2] <- eqdist.etest(C2,sizes=LL,R=B)$p.value
  p2[i,3] <- bd.test(x=A2,y=B2,num.permutations = B,seed=i*12345)$p.value
}
powercomparison2 <- colMeans(p2<alpha)

p3 <- matrix(NA,N,3)
for(i in 1:N){
  A3 <- matrix(rt(L1*p,1),ncol=p);
  B3 <- cbind(rnorm(L2),rnorm(L2,1,1.5));
  C3 <- rbind(A3,B3)
  p3[i,1] <- NNtest(C3,LL,k)$p.value
  p3[i,2] <- eqdist.etest(C3,sizes=LL,R=B)$p.value
  p3[i,3] <- bd.test(x=A3,y=B3,num.permutations = B,seed=i*12345)$p.value
}
powercomparison3 <- colMeans(p3<alpha)

p4 <- matrix(NA,N,3)
L1<-10;L2<-100;L<-L1+L2;LL=c(L1,L2)
for(i in 1:N){
  A4 <- matrix(rnorm(L1*p),ncol=p);
  B4 <- matrix(rnorm(L2*p,1,2),ncol=p);
  C4 <- rbind(A4,B4)
  p4[i,1] <- NNtest(C4,LL,k)$p.value
  p4[i,2] <- eqdist.etest(C4,sizes=LL,R=B)$p.value
  p4[i,3] <- bd.test(x=A4,y=B4,num.permutations = B,seed=i*12345)$p.value
}
powercomparison4 <- colMeans(p4<alpha)
powercomparison <- rbind(powercomparison1,powercomparison2,powercomparison3,powercomparison4)
colnames(powercomparison) <- c('NN','Energy','Ball')
powercomparison

## -----------------------------------------------------------------------------
set.seed(12345)

f <- function(theta,eta,x){   
  stopifnot(theta > 0)
  return( 1/(theta*pi*(1+((x-eta)/theta)^2)) )
}

n <- 2000;
A <- numeric(n);
C <- runif(n);
A[1] <- rnorm(1,0,1);   
i =0;
theta=1;eta=0;

for(j in 2:n){       
  At <- A[j-1]
  B <- rnorm(1,At,1)
  H1 <- f(theta,eta,B)*dnorm(At,B,1)
  H2 <- f(theta,eta,At)*dnorm(B,At,1)
  if(C[j]<=H1/H2)
    A[j] = B
  else{
    A[j] = At
    i = i+1
  }
}

burnin <- 1001 
X <- A[burnin:n]

q1 <- quantile(X,probs = seq(0.1,0.9,0.1))  
q2 <- qt(seq(0.1,0.9,0.1),df=1)   
Q <- rbind(q1,q2)
knitr::kable(head(Q))

## -----------------------------------------------------------------------------
set.seed(12345)
M <- 2*1e3
bi <- 1e3
A <- matrix(0,M,2)
a <- 1 ; b <- 1 ; n <- 10 

A[1,] <- c(0,0)   
for(i in 2:M){    
  y <- A[i-1,2]
  A[i,1] <- rbinom(1,10,y)
  x <- A[i,1]
  A[i,2] <- rbeta(1,x+a,n-x+b)
}

X <- A[(bi+1):M,]    
colMeans(X)  
cov(X)
cor(X)
plot(X,main = "")

## -----------------------------------------------------------------------------

GR <- function(fai){
  fai <- as.matrix(fai)
  n <- ncol(fai)  
  k <- nrow(fai)  
  
  faibar <- rowMeans(fai)
  Bn <- n*var(faibar)
  Wn <- mean(apply(fai, 1, "var"))
  return((Wn*(n-1)/n+(Bn/n))/Wn)
}

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 2000;
m <- 4;

chain93 <- function(){
  A <- numeric(n);
  C <- runif(n);
  A[1] <- 20; 
  i =0;
  theta=1;eta=0;
  
  for(j in 2:n){         
    At <- A[j-1]
    B <- rnorm(1,At,1)
    H1 <- f(theta,eta,B)*dnorm(At,B,1)
    H2 <- f(theta,eta,At)*dnorm(B,At,1)
    if(C[j]<=H1/H2)
      A[j] = B
    else{
      A[j] = At
      i = i+1
    }
  }
  return(A)
}

Z <- matrix(0,m,n)
for(i in 1:m)
  Z[i,] <- chain93()

fai <- t(apply(Z,1,cumsum))
for(i in 1:nrow(fai))
  fai[i,] <- fai[i,]/(1:ncol(fai))
print(GR(fai))

b <- 1000;
rhat <- rep(0,n)
for(j in (b+1):n)
  rhat[j] <- GR(fai[,1:j])
plot(rhat[(b+1:n)], xlab="", ylab="R_hat")
abline(h=1.2)

## -----------------------------------------------------------------------------
set.seed(12345)

M <- 2000;
m <- 4;

chain98 <- function(){
  A <- matrix(0,M,2)
  a <- 1 ; b <- 1 ; n <- 10
  A[1,] <- c(0,0)  
  for(i in 2:M){    
    y <- A[i-1,2]
    A[i,1] <- rbinom(1,10,y)
    x <- A[i,1]
    A[i,2] <- rbeta(1,x+a,n-x+b)
  }
  return(A)
}

Z1 <- matrix(0,m,M)
Z2 <- matrix(0,m,M)
for(i in 1:m){
  Q <- chain98()
  c <- Q[,1]
  d <- Q[,2]
  Z1[i,] <- c
  Z2[i,] <- d
}

fai1 <- t(apply(Z1,1,cumsum))
for(i in 1:nrow(fai1))
  fai1[i,] <- fai1[i,]/(1:ncol(fai1))
print(GR(fai1))

fai2 <- t(apply(Z2,1,cumsum))
for(i in 1:nrow(fai2))
  fai2[i,] <- fai2[i,]/(1:ncol(fai2))
print(GR(fai2))

b <- 1000;
Xrhat <- rep(0,M)
for(j in (b+1):M)
  Xrhat[j] <- GR(fai1[,1:j])
plot(Xrhat[(b+1:M)],xlab="X", ylab="R_hat")
abline(h=1.2)

Yrhat <- rep(0,M)
for(j in (b+1):M)
  Yrhat[j] <- GR(fai2[,1:j])
plot(Yrhat[(b+1:M)], xlab="Y", ylab="R_hat")
abline(h=1.2)

## -----------------------------------------------------------------------------
Yk <- function(k,a){
  yk1 = (-1)^k/factorial(k)/(2^k)
  d = length(a)
  ena = 0           
  for(i in 1:d){
    ena <- ena + a[i]^2
  }
  ena <- sqrt(ena)
  
  yk2 = (ena^(2*k+2))/((2*k+1)*(2*k+2))
  yk3 = exp(lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+d/2+1))
  return(yk1 * yk2 * yk3)
}

## -----------------------------------------------------------------------------
SumofYk <- function(a){
  i=0; sumofyk=0;
  while (abs(Yk(i,a)) > 0.00000001) {    
      sumofyk=sumofyk + Yk(i, a);
      i=i+1;
  } 
  return(sumofyk)
}

## -----------------------------------------------------------------------------
a <- c(1,2)
SumofYk(a)

## -----------------------------------------------------------------------------
aoff <- function(k){
  x <- uniroot( function(a) {pt(sqrt((a^2)*(k-1)/(k-a^2)),df=k-1)-pt(sqrt((a^2)*k/(k+1-a^2)),df=k)}
                ,c(0.00002, sqrt(k)-0.00002))$root
  return(x)
}

N <- c(4:25,100,500,1000)
L <- length(N)
X <- numeric(L)
for(j in 1:L){
  X[j]=aoff(N[j])
}

Y <- cbind(N,X)
colnames(Y) = c("k","a")
Y

## -----------------------------------------------------------------------------
x <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
A <- 3
B <- 10/(sum(x)+3/A)

while (abs(B-A)>0.000000000001){
  A <- B
  B <- 10/(sum(x)+3/A)
}
B

## -----------------------------------------------------------------------------
attach(mtcars)

formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

n1 <- length(formulas)
A1 <- vector("list",n1) 
for(j in seq_along(formulas)){
  A1[[j]] <- lm(formulas[[j]])
}

A2 <- lapply(formulas, lm)
rsq <- function(mod) summary(mod)$r.squared

lapply(A1, rsq)
lapply(A2, rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

n2 <- length(bootstraps)
B1 <- vector("list",length(bootstraps)) 
for(j in seq_along(bootstraps)){
  B1[[j]] <- lm(mpg ~ disp, bootstraps[[j]])
}

B2 <- lapply(bootstraps, lm, formula = mpg ~ disp)
rsq <- function(mod) summary(mod)$r.squared

lapply(B1, rsq)
lapply(B2, rsq)

## -----------------------------------------------------------------------------
pureR <- function(a,b,n,M){
  A <- matrix(0,M,2)
  A[1,] <- c(0,0)   
  for(j in 2:M){     
    y <- A[j-1,2]
    A[j,1] <- rbinom(1,10,y)
    x <- A[j,1]
    A[j,2] <- rbeta(1,x+a,n-x+b)
  }
  return(A)
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('NumericMatrix CppFun(double a, double b, int n, int M) {
            NumericMatrix A(M,2); 
            A(0, 0) = A(0, 1) = 0;
            double x; double y;
            for(int j = 1; j < M; j++)
            {
              y = A(j-1, 1);
              A(j, 0) = rbinom(1, n, y)[0];
              x = A(j, 0);
              A(j, 1) = rbeta(1, x+a, n-x+b)[0];
            }
            return A;
}')

## -----------------------------------------------------------------------------
set.seed(12345)
library(Rcpp)
library(microbenchmark)
a <- 1; b <- 1; n <- 10; M <- 5*1e3;     

Z1 <- pureR(1,1,10,5*1e3)       
Z2 <- CppFun(1,1,10,5*1e3)

par(mfrow = c(1, 2))       
qqplot(Z1[1001:M,1],Z2[1001:M,1],xlab='pure R language_x',ylab='Rcpp function_x')
qqplot(Z1[1001:M,2],Z2[1001:M,2],xlab='pure R language_y',ylab='Rcpp function_y')

## -----------------------------------------------------------------------------
library(microbenchmark)
timecompare <- microbenchmark(R=pureR(1,1,10,5*1e3),Rcpp=CppFun(1,1,10,5*1e3))  
print(summary(timecompare)[, c(1,3,5,6)])

