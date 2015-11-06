library(rstan)

# linear regression with normal errors and partially observed outcomes

N <- 1000
K <- 3
N_mis <- 0.975 * N
N_obs <- N - N_mis

x <- cbind(1, matrix(rnorm(N*K), ncol = K))

beta <- runif(K + 1, 0, 10)

sigma <- abs(rnorm(1))

y <- x %*% beta + rnorm(N, sd = sigma)
y <- as.numeric(y)

y_obs <- y[1:N_obs]
y_mis <- y[(N_obs + 1):N]

x_obs <- x[1:N_obs, ]
x_mis <- x[(N_obs + 1):N, ]

model_string <- "
data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  int<lower=0> K;
  vector[N_obs] y_obs;
  matrix[N_obs, K] x_obs;
  matrix[N_mis, K] x_mis;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
  vector[N_mis] y_mis;
}
model {
  y_obs ~ normal(x_obs * beta, sigma);
  y_mis ~ normal(x_mis * beta, sigma);
}
"

stan_data <- list(N_obs = N_obs, 
                  N_mis = N_mis, 
                  K = K + 1, 
                  y_obs = y_obs, 
                  x_obs = x_obs, 
                  x_mis = x_mis)

str(stan_data)

fit1 <- stan(model_code = model_string, data = stan_data, iter = 5000)

fit1

beta
sigma