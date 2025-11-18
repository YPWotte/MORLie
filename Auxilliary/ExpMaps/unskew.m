function [x] = unskew(A)
%UNSKEW Unskews matrices in R^3x3 or R^4x4, assumed to be in tilde form 
% Idea is that unskew(skew(x)) = x, with skew as defined in skew.m
% Also does unskew(ad(x)) = x.
if length(A(1,:)) == 3 && length(A(:,1)) == 3
    x = [A(3,2); A(1,3); A(2,1)];
elseif length(A(1,:)) == 4 && length(A(:,1)) == 4
    x = unskew(A(1:3,1:3));                     
    x(4:6) = A(1:3,4);
elseif length(A(1,:))==6 && length(A(:,1)) == 6
    x = unskew(A(1:3,1:3));
    x(4:6) = unskew(A(4:6,1:3));
end  
end

