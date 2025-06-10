# LMM Kernel Simulation
#-------------------------------------------#
# Simulation 1 inner product kernel
#-------------------------------------------#
source("./code/REML-EM.R")

variant<-function(x)
{
  length(unique(x))>1
}

inner_product_kernel <- function(n, out=NULL)
{
  p <- 2*n
  allele_freq <- runif(p, 0.001, 0.02)
  
  # generate genotype matrix using Hardy-Weinberg equilibrium
  G <- matrix(0, n, p)
  for(i in 1:p)
  {
    G[,i] <- sample(c(0,1,2), size = n, replace = TRUE, prob = 
                    c((1-allele_freq[i])^2, 2*allele_freq[i]*(1-allele_freq[i]), allele_freq[i]^2))
  }
  vrt <- apply(G, 2, variant)
  G <- as.matrix(G[,vrt])
  G.std <- apply(G, 2, scale)
  
  W <- diag(dbeta(allele_freq[vrt], 1, 25)^2)
  W <- W / max(W)
  K <- 1/p * G.std %*% W %*% t(G.std)
  
  # generate design matrix for fixed effects
  X <- cbind(1, matrix(rnorm(2*n), n, 2))

  # generate response
  a <- mvrnorm(1, mu = rep(0,n), Sigma = 0.6*K)
  e <- rnorm(n, 0, sqrt(2))
  b <- c(1,2,-1)
  y <- X %*% b + a + e

  #vc_list <- LMM_REML(y = y, X = X, beta = b, theta = c(1, 1), K = K)
  vc_list <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K)
  theta <- vc_list$Theta[vc_list$num_iter,]
  
  if(!is.null(out))
            saveRDS(theta, out)
  invisible(theta)
  
}







