function [FA] = Cayley_Hamilton_Function(p, skewA, A)
%CAYLEY_HAMILTON_FUNCTION Summary of this function goes here

pA = p(A);
n = length(pA);
FA = pA(1)*skewA^0;
for i = 2:n
    FA = FA + pA(i)*skewA^(i-1); 
end

end

