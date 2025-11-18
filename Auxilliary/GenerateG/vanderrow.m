function [row] = vanderrow(v,n,m)
%VANDERROW Row of Vandermonde matrix to nth order, for element v
%   Confluent to m-th order, if required
if ~exist('m','var')
    m = 0;
end

if ~exist('n','var')
    n = 1;
end

for i = 0:n-1
    if i >= m
        row(i+1) = prod((i+1-m):i)*v^(i-m);
    end
end

end

