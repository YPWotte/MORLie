function [p] = Cayley_Hamilton_Coefficients(f,A)
%CAYLEY_HAMILTON_COEFFICIENTS returns Cayley-Hamilton Coefficients to turn
%scalar function f into matrix function f(A) = sum_i=0^n-1 p(A)_i A^i.
%   Assumes symbolic f, such that f(x,k) = (d^k/dx^k f)(x)

lambda = eig(A);
n = length(lambda);

if isa(A,'double')
    trunc = 1e15;
    lambda = (round(real(lambda*trunc))+1i*round(imag(lambda*trunc)))/trunc;
end

[GC,GR] = groupcounts(lambda);
nGC = length(GC);

k = 0;
for i = 1:nGC
    for j = 1:GC(i)
        k = k+1;
        M(k,:) = vanderrow(GR(i),n,j-1);
        F(k,1) = f(GR(i),j-1);
    end
end
p = M\F;



end

% matlabFunction(p,"File","matrixF","Vars",{T});