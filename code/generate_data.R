## Generate a population of size 10000 to sample from
## p=m*10000 being the number of genetic variants
## sigma_snp is used to generate the bi-directional genetic effect from a normal distribution
## causal is the causal SNP ratio
## effect is used to specify the nonlinear genetic relation

generate_data <- function(m=2, sigma_snp = 0.2, causal = 1, effect = c("linear", "quad", "trig"))
{
    n <- 10000
    p <- m*n
    allele_freq <- runif(p, 0.05, 0.5)
    
  
    # generate genotype matrix using Hardy-Weinberg equilibrium
    G <- matrix(0, n, p)
    for(i in 1:p)
    {
        G[,i] <- sample(c(0,1,2), size = n, replace = TRUE, prob = 
                    c((1-allele_freq[i])^2, 2*allele_freq[i]*(1-allele_freq[i]), allele_freq[i]^2))
    }
    
    # generate design matrix for fixed effects
    X <- cbind(1, matrix(rnorm(2*n), n, 2))
    b <- c(1,2,-1)

    # generate response
    e <- rnorm(n, 0, sqrt(2))
    causal_loci <- sample.int(p, causal * p)
    if(effect == "linear")
    {
        gno_effect <- G[, causal_loci] %*% rnorm(length(causal_loci), 0, sigma_snp)
    }
  
    if(effect == "quad")
    {
        gno_effect <- (G[, causal_loci] %*% rnorm(length(causal_loci), 0, sigma_snp))^2
    }
  
    if(effect == "trig")
    {
        gno_effect <- cos(G[, causal_loci] %*% rnorm(length(causal_loci), 0, sigma_snp))
    }
    y <- X %*% b + gno_effect + e
    herit <- var(gno_effect) / (2+var(gno_effect))
    
    return_list <- list(heritability = herit, y = y, X = X, genotype_matrix = G)
    return(return_list)
}


## run the following lines of code to generate data and saveRDS
## data_list <- generate_data(sigma_snp = 0.01, causal = 0.4, effect = "quad") 
## saveRDS(data_list, "data_list_sigmaSNP01_causal4_quad.rds")

## data_list <- generate_data(sigma_snp = 0.05, causal = 0.4, effect = "trig") 
## saveRDS(data_list, "data_list_sigmaSNP05_causal4_trig.rds")


## data_list <- generate_data(sigma_snp = 0.01, causal = 0.6, effect = "quad") 
## saveRDS(data_list, "data_list_sigmaSNP01_causal6_quad.rds")

## data_list <- generate_data(sigma_snp = 0.05, causal = 0.6, effect = "trig") 
## saveRDS(data_list, "data_list_sigmaSNP05_causal6_trig.rds")

## data_list <- generate_data(sigma_snp = 0.01, causal = 0.8, effect = "quad") 
## saveRDS(data_list, "data_list_sigmaSNP01_causal8_quad.rds")

## data_list <- generate_data(sigma_snp = 0.05, causal = 0.8, effect = "trig") 
## saveRDS(data_list, "data_list_sigmaSNP05_causal8_trig.rds")


