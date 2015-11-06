library(rstan)

# multivariate normal model

model_string <- "data {
  int<lower=0> N;
  vector[2] y[N];
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
}
model {
  y ~ multi_normal(mu, Sigma);
}"

y <- MASS::mvrnorm(1000, c(0, 0), matrix(c(1, 0.5, 0.5, 2), nrow = 2))

plot(y[, 1], y[, 2])

mydata <- list(y, 
               N = nrow(y), 
               var1 = var(y[, 1]), 
               var2 = var(y[, 2]))

fit1 <- stan(model_code = model_string, data = mydata)

fit1