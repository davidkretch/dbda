# find the sample size necessary to reject null with probability p

library(rstan)

model_code <- "
data {
int<lower=0> N;
real y[N];
}
parameters {
real mu;
real<lower=0> sigma;
}
model {
y ~ normal(mu, sigma);
}
"

# build the model; do not estimate
fit <- stan(model_code = model_code, 
            data = list(N = 2, y = c(0, 1)), 
            chains = 1, iter = 1)

hdi <- function(sampleVec, credMass=0.95) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

sample_sizes <- c(50, 100, 200)
iter <- 1000
sim_results <- matrix(nrow = length(sample_sizes), ncol = iter)
rownames(sim_results) <- paste0("n", sample_sizes)

overlap <- function(x, y) {x[1] <= y[2] & y[1] <= x[2]}

# test each sample size 1000 times
for (i in seq_along(sample_sizes)) {
  for (j in seq(iter)) {
    sample_size <- sample_sizes[i]
    
    # 1. generate representative value
    mu <- rnorm(1)
    
    # 2. generate random sample
    y <- rnorm(sample_size, mean = mu)
    data <- list(N = sample_size, y = y)
    
    # 3. estimate posterior
    posterior <- stan(data = data, fit = fit, chains = 1)
    
    # 4. does HDI exclude ROPE?
    hdi_mu <- hdi(extract(posterior)$mu)
    success <- !overlap(c(-0.02, 0.02), hdi_mu)
    sim_results[i, j] <- success
  }
}

# power: probability of rejecting the null, by sample size
apply(sim_results, 1, mean)