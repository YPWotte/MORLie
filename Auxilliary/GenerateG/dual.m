function [dEi] = dual(Ei)
%DUAL Returns dual of basis Ei, assumed to be a subset of m x m

n = size(Ei,1);
[m,~] = size(Ei{1});
N = 10; % precision
Ei_flat = zeros(m*m);
dEi = cell(n,1);
for i = 1:n
    Ei_flat(i,:) = reshape(Ei{i},[],1);
end
Ei_flat = [Ei_flat(1:n,:); reduce_base(Ei_flat.').'];

dEi_flat = round(inv(Ei_flat).',N);
for i = 1:n
    dEi{i} = reshape(dEi_flat(i,:),m,m);
end

end

