function X = simpleTruncatedNorm(sigma, a, b, n, mu)
if nargin < 5
  mu = 0;
end

PhiA = normcdf(a, mu, sigma);
PhiB = normcdf(b, mu, sigma);
X = rand(1,n);
X = PhiA + X*(PhiB-PhiA);
X = norminv(X, mu, sigma);

