library(rstan)

set.seed(1)
y <- rbinom(100, size = 1, prob = 0.75)

omega1 <- .25
kappa1 <- 12
omega2 <- .75
kappa2 <- 12

a1 <- omega1*(kappa1-2)+1
b1 <- (1-omega1)*(kappa1-2)+1
a2 <- omega2*(kappa2-2)+1
b2 <- (1-omega2)*(kappa2-2)+1

# plot the priors
curve(dbeta(x, shape1 = a1, shape2 = b1), 0, 1)
curve(dbeta(x, shape1 = a2, shape2 = b2), 0, 1)

model_string <- "
data {
  int<lower=0> N;
  int<lower=0,upper=1> y[N];
}
parameters {
  real<lower=0,upper=1> m;
  real<lower=0,upper=1> theta1;
  real<lower=0,upper=1> theta2;
}
transformed parameters {
  real<lower=0,upper=1> theta;
  theta <- if_else(m < 0.5, theta1, theta2);
}
model {
  m ~ uniform(0, 1);
  theta1 ~ beta(3.5, 8.5);
  theta2 ~ beta(8.5, 3.5);
  y ~ bernoulli(theta);
}
"

stan_data <- list(N = length(y), y = y)

fit <- stan(model_code = model_string, data = stan_data, iter = 1e4)

fit

plot(fit)

prop.table(table(round(fit@sim$samples[[1]]$m)+1))
