function [inner] = duality_product(Ei,dEj)
%DUALITY_PRODUCT Duality product of Ei and dual element dEj
%   Expects Ei and dEj to be m x m matrices.
Ei_flat = reshape(Ei,[],1);
dEi_flat = reshape(dEj,[],1);
inner = dEi_flat.'*Ei_flat;
end

