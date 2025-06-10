source("./code/REML-EM.R")


narrow_heritability <- function(n, tag=1, c=1, power = 2, sigma = NULL, out=NULL)
{
  if(tag == 1)
  {
      data_list <- readRDS("./code/data_list_sigmaSNP01_causal8_quad.rds")
  }
  
  if(tag == 2)
  {
      data_list <- readRDS("./code/data_list_sigmaSNP05_causal8_trig.rds")
  }
  
  herit <- data_list$heritability
  
  ## sample n subjects from the population
  id <- sample(10000, n)
  y <- data_list$y[id]
  X <- data_list$X[id,]
  G <- data_list$genotype_matrix[id,]
  G.std <- apply(G, 2, scale)
  p <- ncol(G.std)
  
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