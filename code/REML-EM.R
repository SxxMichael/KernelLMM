# This algorithm provides REML estimate of an LMM of the form using EM algorithm
# y=Xb+a+e, a~N(0,sigma_a^2 * K), e~N(0, sigma_e^2 * I)
# y is the vector of response
# X is the design matrix for fixed effect
# beta is the initial value for fixed effect
# theta = (sigma_a^2, sigma_e^2) is the intial value for variance component 
# K is the covariance matrix of the random effect
library(Matrix)
library(MASS)

# trace function
tr <- function(A)
{
  return(sum(diag(A)))
}

LMM_REML <- function(y, X, theta, K, num_iter = 1e3, tol = 1e-7)
{
  n <- length(y)
  #q <- rankMatrix(X)
  
  id_mat <- diag(rep(1,n))
  
  
  # matrix for parameter updates in EM algorithm
  Theta <- matrix(0, num_iter, 2)
  Theta[1,] <- theta
  
  for(i in 2:num_iter)
  {
    sigma_a.i <- Theta[i-1,1]
    sigma_e.i <- Theta[i-1,2]
    V.i <- sigma_a.i * K + sigma_e.i * id_mat
    V_inv.i <- solve(V.i)
    P.i <- V_inv.i - V_inv.i %*% X %*% solve(t(X) %*% V_inv.i %*% X) %*% t(X) %*% V_inv.i
    
    # update parameters
    sigma_a.i1 <- sigma_a.i + 
                  (1/n) * sigma_a.i^2 * t(y) %*% P.i %*% K %*% P.i %*% y -
                  (1/n) * sigma_a.i^2 * tr(P.i %*% K)
    Theta[i,1] <- sigma_a.i1
    
    sigma_e.i1 <- sigma_e.i +
                  (1/n) * sigma_e.i^2 * t(y) %*% P.i %*% P.i %*% y -
                  (1/n) * sigma_e.i^2 * tr(P.i)
    Theta[i,2] <- sigma_e.i1
                  
    
   # gamma.i   <- sigma_a.i/sigma_e.i
    
    # error contrast matrix
    #qrX <- qr(X)
    #y.tilde <- qr.qty(qrX,y)[(rankMatrix(X)+1):n]
    #A.t <-t(qr.Q(qrX,complete=TRUE)[,(rankMatrix(X)+1):n])
    #A <- t(A.t)
    
    # marginal covariance matrix for y
    #V_gamma <- A.t %*% K %*% A
    #Sigma_gamma <- gamma.i * V_gamma + id_mat
    #Sigma_gamma.inv <- solve(Sigma_gamma)
    
    # update parameters
    #sigma_a.i1 <- sigma_a.i + 
    #              1/(n-q) * gamma.i^2 * t(y.tilde) %*% Sigma_gamma.inv %*% V_gamma %*% Sigma_gamma.inv %*% y.tilde -
    #              1/(n-q) * sigma_a.i^2/sigma_e.i * tr(Sigma_gamma.inv %*% V_gamma)
    #Theta[i,1] <- sigma_a.i1
    
    #sigma_e.i1 <- 1/(n-q) * (sigma_a.i * tr(V_gamma) - sigma_a.i^2/sigma_e.i * tr(V_gamma %*% Sigma_gamma.inv %*% V_gamma)
    #                         + crossprod(t(id_mat - gamma.i * Sigma_gamma.inv %*% V_gamma) %*% y.tilde))
    #Theta[i,2] <- sigma_e.i1
    
    #cat(Theta[i,], "\n")
    
    if(crossprod(Theta[i,]-Theta[i-1,]) < tol)
    {
      break
    }
    
  }
  
  return_list <- list(Theta = Theta, num_iter = i)
  return(return_list)
}