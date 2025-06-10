source("./code/REML-EM.R")


narrow_heritability <- function(n, m=4, sigma_snp = 0.1, causal = 0.5, c=1, power = 2, sigma = NULL, out=NULL)
{
  p <- m*n
  allele_freq <- runif(p, 0.05, 0.5)
  
  # generate genotype matrix using Hardy-Weinberg equilibrium
  G <- matrix(0, n, p)
  for(i in 1:p)
  {
    G[,i] <- sample(c(0,1,2), size = n, replace = TRUE, prob = 
                    c((1-allele_freq[i])^2, 2*allele_freq[i]*(1-allele_freq[i]), allele_freq[i]^2))
  }
  G.std <- apply(G, 2, scale)
  
  K_prod <- (1/p) * tcrossprod(G.std)
  K_poly <- (c + K_prod)^2
  
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
  

  
  # generate design matrix for fixed effects
  X <- cbind(1, matrix(rnorm(2*n), n, 2))
  b <- c(1,2,-1)

  # generate response
  e <- rnorm(n, 0, sqrt(2))
  causal_loci <- sample.int(p, causal * p)
  gno_effect <- G.std[, causal_loci] %*% rnorm(length(causal_loci), 0, sigma_snp)
  y <- X %*% b + gno_effect + e
  herit <- var(gno_effect) / var(y)

  vc_list_prod <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_prod)
  vc_list_poly <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_poly)
  vc_list_gaussian <- LMM_REML(y = y, X = X, theta = c(1, 1), K = K_gaussian)
  
  theta_prod <- vc_list_prod$Theta[vc_list_prod$num_iter,]
  theta_poly <- vc_list_poly$Theta[vc_list_poly$num_iter,]
  theta_gaussian <- vc_list_gaussian$Theta[vc_list_gaussian$num_iter,]
  
  herit_prod <- theta_prod[1]/sum(theta_prod)
  herit_poly <- theta_poly[1]/sum(theta_poly)
  herit_gaussian <- theta_gaussian[1]/sum(theta_gaussian)
  
  theta <- c(herit, herit_prod, herit_poly, herit_gaussian)
  print(theta)
  
  if(!is.null(out))
            saveRDS(theta, out)
  invisible(theta)

}  