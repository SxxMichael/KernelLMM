source("./code/REML-EM.R")


mis_kernel <- function(n, tag=1, c=1, power = 2, sigma = NULL, out=NULL)
{
  p <- 2*n
  allele_freq <- runif(p, 0.05, 0.5)
  
  # generate genotype matrix using Hardy-Weinberg equilibrium
  G <- matrix(0, n, p)
  for(i in 1:p)
  {
    G[,i] <- sample(c(0,1,2), size = n, replace = TRUE, prob = 
                    c((1-allele_freq[i])^2, 2*allele_freq[i]*(1-allele_freq[i]), allele_freq[i]^2))
  }
  G.std <- apply(G, 2, scale)
  
  dist_vec <- c(dist(G.std))
  
  if(is.null(sigma))
  {
    sigma_sq <- var(dist_vec)
  }else
  {
    sigma_sq <- sigma^2
  }
  K_temp <- as.matrix(dist(G.std)^2)/p
  
  K_gaussian <- exp(-K_temp/(2*sigma_sq))
  K_prod <- tcrossprod(G.std)/p
  K_poly <- (c+K_prod)^power
  
  # generate design matrix for fixed effects
  X <- cbind(1, matrix(rnorm(2*n), n, 2))

  # generate response
  # tag == 1, truth = prod
  # tag == 2, truth = polynomial
  # tag == 3, truth = gaussian
  e <- rnorm(n, 0, sqrt(2))
  b <- c(1,2,-1)
  if(tag == 1)
  {
    a <- mvrnorm(1, mu = rep(0,n), Sigma = 0.6*K_prod)
    y <- X %*% b + a + e

    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_poly)
    theta_poly <- vc_list$Theta[vc_list$num_iter,]
    
    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_gaussian)
    theta_gaussian <- vc_list$Theta[vc_list$num_iter,]
    
    theta <- c(theta_poly, theta_gaussian)
  
  }
  
  if(tag == 2)
  {
    a <- mvrnorm(1, mu = rep(0,n), Sigma = 0.6*K_poly)
    y <- X %*% b + a + e
    
    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_prod)
    theta_prod <- vc_list$Theta[vc_list$num_iter,]
    
    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_gaussian)
    theta_gaussian <- vc_list$Theta[vc_list$num_iter,]
    
    theta <- c(theta_prod, theta_gaussian)
  }
  
  if(tag == 3)
  {
    a <- mvrnorm(1, mu = rep(0,n), Sigma = 0.6*K_gaussian)
    y <- X %*% b + a + e
    
    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_poly)
    theta_poly <- vc_list$Theta[vc_list$num_iter,]
    
    vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_prod)
    theta_prod <- vc_list$Theta[vc_list$num_iter,]
    
    theta <- c(theta_prod, theta_poly)
  }
  
  
  
  
  if(!is.null(out))
            saveRDS(theta, out)
  invisible(theta)

}  