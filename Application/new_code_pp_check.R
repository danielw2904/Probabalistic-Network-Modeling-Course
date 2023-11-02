library(data.table)
library(dplyr)
library(ggraph)
library(igraph)
library(countrycode)
library(tidyr)
library(stringr)
library(nimble)
library(bayesplot)

el <- fread("~/Projects/Probabalistic-Network-Modeling-Course/Application/spotify_connections.csv.gz")
el[, streams_a_s := total_streams_a / max(total_streams_a)] # Normalized Streams

el$weight_b <- ifelse(el$weight >= 671/2, 1, 0)

# Remove global region and rename "gb" to "uk"
el <- el %>%
  filter(region_a != "global" & region_b != "global") %>%
  mutate(region_a = str_replace(region_a, "gb", "uk"),
         region_b = str_replace(region_b, "gb", "uk"))

# Use codelist dataframe from countrycode package 
# Remove "." from all names to match those in our dataset
codelist$region <- gsub("\\.", "", codelist$cctld)

# Grab relevant sublist from codelist 
sub_codelist <- codelist[, c("region", "iso3c", "continent", "country.name.en")]

# This creates the node properties
node_props <- el %>% 
  
  # We only need one row per node, so we group by and summarise 
  group_by(region_a) %>%
  summarise(total_streams = min(total_streams_a)) %>%
  
  # Some vertices only show up in region_a and some only in region_b, so we perform a full join to have both
  full_join(data.frame(region = sort(unique(el$region_b))), ., c("region" = "region_a")) %>%
  
  # Merge in the node properties from the codelist
  left_join(sub_codelist, by = "region") %>%
  
  # Replace the NA in argentina with an average (it never appears in the first column, so it has no total_streams_a)
  mutate(total_streams = replace_na(total_streams, mean(total_streams, na.rm = TRUE))) %>%
  
  # Scale down total streams
  mutate(total_streams = total_streams/sd(total_streams))

g <- graph_from_data_frame(el[el$weight_b == 1,], vertices = node_props, directed = FALSE)

adj_mat <- as_adjacency_matrix(g)
A <- as.matrix(adj_mat)


K <- 10
N <- nrow(A)
a_o <- .5
b_o <- 1
a_c <- 1
b_c <- .5

alpha <- 0.01#1/K

# Define the constants
sbmConsts <- list(N = N, K= K,
                  a_o = a_o, b_o = b_o,
                  a_c = a_c, b_c = b_c,
                  alpha =  rep(alpha, K))

# Define the data
sbmData <- list(Y = A)

theta_init <- matrix(.2, K, K)
diag(theta_init) <- 0.6

sbmInits <- list(lambda = rep(1/K, K), # block assignment probs
                 theta = theta_init, # edge probs, better init
                 Pi = matrix(mean(A[lower.tri(A)]), N, N), # edge probs
                 Z = sample(1:K, N, TRUE), # group indicators
                 mat.logL = matrix(0,N,N)) #####ODIFY

# Define model with BUGS code
sbmCode <- nimbleCode({
  
  # lambda - latent group assignment probabilities (vector of length K)
  lambda[1:K] ~ ddirch(alpha[]) 
  
  # Z - latent group indicator (binary vector of length K, summing to 1)
  for(i in 1:N){
    Z[i] ~ dcat(prob = lambda[])
  }
  # theta - symmetric matrix of within and between group edge probabilities
  for (i in 1:K){
    theta[i,i] ~ dbeta(shape1 = a_c, shape2 = b_c) # Within block connections
    for (j in 1:(i-1)){ 
      theta[i,j] ~ dbeta(shape1 = a_o, shape2 = b_o) # Between block connections 
      theta[j,i] <- theta[i,j] # symmetric matrix
    }
  }
  
  # Pi - node to node edge probabilities (based on group membership)
  for (i in 2:N){
    for (j in 1:(i-1)){ # Self-edges not allowed
      Pi[i,j] <- myCalculation(theta[,], Z[i], Z[j]) # Workaround because nimble does not allow indexing by latent variables
      Y[i,j] ~ dbin(size = 1, prob = Pi[i,j])
    }
  }
  
  ##### MODIFY
  # Compute logL
  for (i in 2:N){
    for (j in 1:(i-1)){ # Self-edges not allowed
      mat.logL[i,j] <- log((Pi[i,j]^Y[i,j])*((1 - Pi[i,j])^(1-Y[i,j])))
      #mat.logL[j,i] <- 0
    }}
  
  logL <- sum(mat.logL[1:N, 1:N])/2 # diag is zero so this works
  
})

## User-defined functions: written in NIMBLE
myCalculation <- nimbleFunction(
  run = function(grid = double(2), index1 = double(0), index2 = double(0)) {  ## index could be int() but model variables are represented as double anyway
    return(grid[index1, index2])
    returnType(double(0))
  })

niter <- 500
nburn <- 100
nthin <- 1
nchains <- 2
mcmc.out <- nimbleMCMC(code = sbmCode, constants = sbmConsts,
                       data = sbmData, inits = sbmInits,
                       nchains = nchains, niter = niter, nburnin = nburn, thin = nthin,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('lambda', 'theta', 'Z', "logL", "Pi")) ### MODIFY

mcmc.out$samples <- lapply(mcmc.out$samples, \(x){ x[is.na(x)] <- 0; return(x)}) ### MODIFY
mcmc_trace(mcmc.out$samples, regex_pars = "lambda") ### MODIFY



################### Posterior predictive checks
samples <- do.call(rbind, mcmc.out$samples)

all.samples <- do.call(rbind, mcmc.out$samples) # as matrix
p.samples <- all.samples[, grepl("^Pi", colnames(all.samples))]
p.post <- matrix(data = colMeans(p.samples), byrow = FALSE, nrow = N, ncol = N)
obsA <-  A[upper.tri(A)]
df <- data.frame(pred = p.post[lower.tri(p.post)], data = obsA)
pROC::plot.roc(df$data ~ df$pred, percent = TRUE, print.auc = TRUE, main = "ROC Curve - Baseline GLM")


A.r <- matrix(NA, nrow = N, ncol = N)
postD <- matrix(NA, nrow = nrow(p.samples), ncol = 1)
allD <- matrix(NA, nrow = nrow(p.samples), ncol = N)
for(r in 1:nrow(p.samples)){
  
  # Sample adjacency matrix
  p.post.r <- matrix(data = p.samples[r,], byrow = FALSE, nrow = N, ncol = N)
  diag(p.post.r) <- 0
  A.r[upper.tri(A.r)] <- rbinom(N*(N - 1)/2, 1, p.post.r[lower.tri(p.post.r)])
  A.r[lower.tri(A.r)] <- t(A.r)[lower.tri(A.r)]
  diag(A.r) <- 0
  
  # Sample network density
  postD[r,] <- sum(A.r[lower.tri(A.r)])/(N*(N -1)/2)
  allD[r,] <- rowSums(A.r)
}

D.obs <- sum(A[lower.tri(A)])/(N*(N -1)/2)
hist(postD)
abline(v = D.obs)

all.D.obs <- rowSums(A)
pp_check(all.D.obs, allD[2000:2050, ], ppc_dens_overlay)
