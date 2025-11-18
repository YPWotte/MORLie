function [Eb] = bracket_all(Ei)
%BRACKET_ALL Brackets all vectors in Ei and returns resulting vectors
%   Expects Ei to be an (n,1) cell containing m x m matrices, brackets
%   with bracket on gl(m).
[n,~] = size(Ei);
b = n*(n-1)/2;
Eb = cell(b,1);
k = 0;
for i = 1:n
    for j = i+1:n
        k = k+1;
        Eb{k} = Ei{i}*Ei{j} - Ei{j}*Ei{i};
    end
end

end

