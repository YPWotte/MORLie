function [CA] = adgl(A,C)
%ADGL Returns adjoint representation of twist A, using structure
%coefficients C
%   Expects A to be a column, and C^k_ij collected in C(k,i,j)

if isa(A, "double")
    CA = tensorprod(C,A,2,1);
else
    n = length(A);
    CA = C(:,1,:)*A(1);
    for i = 2:n
        CA = CA + C(:,i,:)*A(i);
    end
    CA = squeeze(CA);
end

end

