function sigma = findSigma(mu, a, b, p)
sigma = max(abs([a b]-mu))*.682/10 * p
tol = 1e-6;
c1 = mu-max([abs(mu-a) abs(mu-b)])
c2 = mu+max([abs(mu-a) abs(mu-b)])

%sigma = fzero(@(x) normcdf(c2,mu,x)-normcdf(c1,mu,x)-p, [0 10e6]);
sigma = fzero(@(x) normcdf(c2,mu,x)-normcdf(c1,mu,x)-p, [0 10e6]);



