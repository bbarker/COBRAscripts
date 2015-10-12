function rvec = randNoResample(rmax, rlen, excludes)
%This will work best with smaller ranges.

if rmax < rlen
  disp('Warning: sample pool smaller than requested sample size.')
  rvec = 1:rmax;
  return
end

if nargin < 3
  excludes = [];
end


samplevec = 1:rmax;

rvec = randi(rmax,1,rlen);
rvec = setdiff(rvec,excludes);

lrvec = length(rvec);
while lrvec < rlen
  rvecD = randi(rmax,1,rlen - lrvec);
  rvec = setdiff(union(rvec, rvecD),excludes);
  lrvec = length(rvec);
end