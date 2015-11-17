library(rstan)

set.seed(1)

N <- 1000
K <- 6
L <- 2

X <- 10 * matrix(runif(N * (K - 1)), ncol = (K - 1))
X <- cbind(1, X)

Beta <- 10 * matrix(rnorm(K * L), ncol = L)

mat <- matrix(runif(L * L), ncol = L)
Sigma <- mat %*% t(mat)

Y <- X %*% Beta + MASS::mvrnorm(N, mu = rep(0, L), Sigma = Sigma)

stan_data <- list(N = N, 
                  K = K, 
                  L = L, 
                  X = X, 
                  Y = Y)

model_string <- "
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> L;
  matrix[N,K] X;
  matrix[N,L] Y;
}
parameters {
  matrix[K,L] Beta;
  cov_matrix[L] Sigma;
}
transformed parameters {
  matrix[N,L] Yhat;
  Yhat <- X * Beta;
}
model {
  for (i in 1:N)
    Y[i] ~ multi_normal(Yhat[i], Sigma);
}
"

fit <- stan(model_code = model_string, data = stan_data)

params <- extract(fit)

Beta_params <- params$Beta
Sigma_params <- params$Sigma

fit_summary <- print(fit)[1:20, ]

dimnames(Beta_params)
Beta

Sigma_params
Sigma
