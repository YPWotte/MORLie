function [C] = structure_coeffs(Ei)
%STRUCTURE_COEFFS Returns structure coefficients of Lie algebra Ei
%   Expects Ei to be an (n,1) cell containing m x m matrices forming a full
%   basis for a Lie algebra. Structure coefficients w.r.t. bracket on gl(m). 

[n,~] = size(Ei);
dEi = dual(Ei);
C = zeros(n,n,n);

for i = 1:n
    for j = i+1:n
        Eij = Ei{i}*Ei{j} - Ei{j}*Ei{i};
        for k = 1:n
            cijk = duality_product(Eij,dEi{k});
            C(k,i,j) = cijk;
            C(k,j,i) = -cijk;
        end
    end
end

end

