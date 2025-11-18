function [skewA] = skewgl(A,Ei)
%SKEWGL Returns skew-form of column A, optionally w.r.t. basis Ei
%   Without Ei, A is padded with zeros to return  a square matrix of 
%   minimal dimension.  
[n,~] = size(A);
if ~exist('Ei','var')
    m = ceil(sqrt(n));
    A = [A; zeros(m^2-n,1)];
    skewA = reshape(A,m,m);
else
    skewA = A(1)*Ei{1};
    for i = 2:n
        skewA = skewA + A(i)*Ei{i};
    end
end

end

