library(rstan)

y1 <- c(70,80,79,83,77,75,84,78,75,75,78,82,74,81,72,70,75,72,76,77)
y2 <- c(56,80,63,62,67,71,68,76,79,67,76,74,67,70,62,65,72,72,69,71)

model_string <- "
data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector[n1] y1;
  vector[n2] y2;
}
parameters {
  real mu1;
  real mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}
transformed parameters {
  real<lower=0> diff;
  diff <- fabs(mu1 - mu2);
}
model {
  y1 ~ normal(mu1, sigma1);
  y2 ~ normal(mu2, sigma2);
}
"

stan_data <- list(n1 = length(y1), 
                  n2 = length(y2), 
                  y1 = y1, 
                  y2 = y2)

fit <- stan(model_code = model_string, data = stan_data)

fit

library(ggplot2)

ggplot() + 
  geom_density(aes(y1)) + 
  geom_vline(xintercept = mean(y1)) + 
  geom_density(aes(y2)) + 
  geom_vline(xintercept = mean(y2))

diff_samples <- extract(fit)$diff

ggplot() + 
  geom_histogram(aes(diff_samples))
