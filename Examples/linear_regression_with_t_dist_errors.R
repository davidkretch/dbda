library(rstan)

# linear regression with t-distributed errors
# gamma(2, 0.1) prior on t-distribution DOF parameter, as suggested by 
# http://andrewgelman.com/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/

N <- 1000
K <- 3

x <- cbind(1, matrix(rnorm(N*K), ncol = K))

beta <- runif(K + 1, 0, 10)

sigma <- abs(rnorm(1))

y <- x %*% beta + rnorm(N, sd = sigma)
y <- as.numeric(y)

# plot(x[, 2], y)

model_string <- "
data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] x;
  vector[N] y;
}
parameters {
  vector[K] beta;
  real<lower=0> nu;
  real<lower=0> sigma;
}
model {
  nu ~ gamma(2, 0.1);
  y ~ student_t(nu, x * beta, sigma);
}
"

stan_data <- list(N = nrow(x), K = ncol(x), x = x, y = y)

str(stan_data)

fit1 <- stan(model_code = model_string, data = stan_data)

fit1
beta
sigma