data {
  int<lower=0> n1;
  int<lower=0> n2;
  int<lower=0> N;
  int<lower=1,upper=n1> blk1[N];
  int<lower=1,upper=n2> blk2[N];
  int<lower=0> y[N];
}
transformed data {
  vector[n1] alpha_a = rep_vector(5, n1);
  vector[n2] alpha_b = rep_vector(5, n2);
}
parameters {
  simplex[n1] a_norm;
  simplex[n2] b_norm;
  real<lower=0> mu;
}
model {
  vector[n1] scaled_a_norm = mu * n1 * n2 * a_norm;
  matrix[n1, n2] gamma;
  vector[N] lin_pred;
  for (j in 1:n2)
    gamma[ , j] = scaled_a_norm * b_norm[j];
  for (n in 1:N)
    lin_pred[n] = gamma[blk1[n], blk2[n]];

  a_norm ~ dirichlet(alpha_a);
  b_norm ~ dirichlet(alpha_b);

  mu ~ gamma(2, 0.1);
  y ~ poisson(lin_pred);
}
generated quantities {
  vector[n1] a = n1 * a_norm;
  vector[n2] b = n2 * b_norm;
}
