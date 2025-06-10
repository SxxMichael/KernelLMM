source("./code/REML-EM.R")

newton_raphson <- function(y, X, beta, theta, K, num_iter = 1e3, tol = 1e-7)
{
    n <- length(y)
    q <- rankMatrix(X)
    id_mat <- diag(rep(1,n))
    Theta <- matrix(0, num_iter, 2)
    Theta[1,] <- theta
    
    for(i in 1:(num_iter-1))
    {
        gamma <- theta[1]/theta[2]
        V_mat <- id_mat + gamma * K
        V_inv <- solve(V_mat)
        P_mat <- V_inv - V_inv %*% X %*% solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
        PKP <- P_mat %*% K %*% P_mat
        
        D2 <- -(n-q)/(2*theta[2]) + 1/(2*theta[2]^2) * t(y) %*% P_mat %*% y
        D1 <- -(1/2)*tr(P_mat %*% K) + 1/(2*theta[2]) * t(y) %*% PKP %*% y
        D_vec <- c(D1, D2)
        
        I11 <- (1/2) * tr(PKP %*% K)
        I12 <- tr(PKP %*% V_mat)/(2*theta[2])
        I22 <- (n-q)/(2*theta[2]^2)
        Fisher_info <- matrix(c(I11, I12, I12, I22), 2, 2)
        
        Theta[i+1,] <- Theta[i,] + solve(Fisher_info, D_vec)
        theta <- Theta[i+1,]
        cat(theta, "\n")
        
        if(crossprod(Theta[i+1,]-Theta[i,]) < tol)
        {
          break
        }
        
    }
    
    return_list <- list(Theta = Theta, num_iter = i)
    return(return_list)
}


vc_est <- function(y, X, beta, theta, K, num_iter = 1e3, tol = 1e-3)
{
    em_list <- LMM_REML(y=y, X=X, beta=beta, theta=theta, K=K, num_iter=num_iter, tol=tol)
    em_i <- em_list$num_iter
    em_est <- em_list$Theta[em_i,]
    
    nr_list <- newton_raphson(y=y, X=X, beta=beta, theta=em_est, K=K)
    return(nr_list)
}
