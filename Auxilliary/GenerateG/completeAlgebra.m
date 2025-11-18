function [Ei_Full] = completeAlgebra(Ei)
%COMPLETEALGEBRA Returns Lie algebra with germ Ei
%   Expects Ei to be an (n,1) cell containing m x m matrices, brackets
%   with bracket on gl(m).

[n,~] = size(Ei);
[m,~] = size(Ei{1});
converged = false;
n_old = n;

while ~converged
    Eb = bracket_all(Ei);
    E_big = [Ei;Eb];
    [n_big,~] = size(E_big);
    Ei_flat = zeros(n_big,m*m);
    for i = 1:n_big
        Ei_flat(i,:) = reshape(E_big{i},[],1);
    end
    Ei_red = reduce_base(reduce_base(Ei_flat.')).';
    [n_new,~] = size(Ei_red);
    Ei_new = cell(n_new,1);
    for i = 1:n_new
        Ei_new{i} = reshape(Ei_red(i,:),m,m);
    end
    if n_new == n_old
        converged = true;
        Ei_Full = Ei_new;
    else
        Ei = Ei_new;
        n_old = n_new;
    end
end

end

