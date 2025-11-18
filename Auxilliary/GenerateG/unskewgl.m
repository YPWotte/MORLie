function [A] = unskewgl(skewA,dEi)
%UNSKEWGL Undoes skew-form, returning a column A
%   Optionally w.r.t. supplied dual basis dEi. 
%   Without Ei, skewA is returned as a column. Trailing zeros are kept.

A = reshape(skewA,[],1);

if exist('dEi','var')
    [n,~] = size(dEi);
    [m,~] = size(dEi{1});
    dEi_flat = zeros(n,m^2);
    for i = 1:n
        dEi_flat(i,:) = reshape(dEi{i},1,[]);
    end
    A = dEi_flat*A;
end

