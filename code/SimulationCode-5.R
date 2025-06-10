source("./code/REML-EM.R")


inner_product_kernel <- function(n, c=1, power = 2, out=NULL)
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
  K <- (c + (1/p) * tcrossprod(G.std))^2
  
  # generate design matrix for fixed effects
  X <- cbind(1, matrix(rnorm(2*n), n, 2))

  # generate response
  a <- mvrnorm(1, mu = rep(0,n), Sigma = 0.6*K)
  e <- rnorm(n, 0, sqrt(2))
  b <- c(1,2,-1)
  y <- X %*% b + a + e

  vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K)
  theta <- vc_list$Theta[vc_list$num_iter,]
  
  if(!is.null(out))
            saveRDS(theta, out)
  invisible(theta)
  
}

