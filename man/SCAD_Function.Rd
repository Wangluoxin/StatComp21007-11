% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCAD_Function.R
\name{SCAD_Function}
\alias{SCAD_Function}
\title{SCAD penalty iterative algorithm}
\usage{
SCAD_Function(X, y)
}
\arguments{
\item{X}{Design Matrix}

\item{y}{Dependent variable}
}
\value{
Regression coefficient and intercept.
}
\description{
This method is used for variable selection and parameter determination for classical linear model.
}
\examples{
\dontrun{
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
beta<-c(3,1.5,0,0,2,0,0,0)
y<- X \%*\% beta + rnorm(n = 60, 0, 1)
SCAD_Function(X,y)
}
}
