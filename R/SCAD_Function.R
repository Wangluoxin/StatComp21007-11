#' @title SCAD penalty iterative algorithm
#' @description This method is used for variable selection and parameter determination for classical linear model.
#' @param X Design Matrix
#' @param y Dependent variable
#' @return Regression coefficient and intercept.
#' @examples
#' \dontrun{
#' Simu_Multi_Norm<-function(x_len, sd = 1, pho = 0.5){
#'   V <- matrix(data = NA, nrow = x_len, ncol = x_len)
#'   for(i in 1:x_len){ 
#'     for(j in 1:x_len){ 
#'       V[i,j] <- pho^abs(i-j)
#'     }
#'   }
#'   V<-(sd^2) * V
#'   return(V)
#' }
#' set.seed(123)
#' X<-mvrnorm(n = 60, mu = rep(0,8), Simu_Multi_Norm(x_len = 8,sd  = 1, pho = 0.5))
#' beta<-c(3,1.5,0,0,2,0,0,0)
#' y<- X %*% beta + rnorm(n = 60, 0, 1)
#' SCAD_Function(X,y)
#' }
#' @export
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

