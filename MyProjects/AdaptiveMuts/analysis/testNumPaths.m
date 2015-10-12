function pathsum = testNumPaths(n)
%Counts the number of edges, not paths
pathsum = 0;
for i = 1:n
 pathsum = pathsum + (n-i+1)*nchoosek(n,i-1);
end
