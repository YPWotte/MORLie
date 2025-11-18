function [Hinv] = inv_se3(H)
%INVH Quick inverse of H in SE(3)

Rt = H(1:3,1:3).';
p = H(1:3,4);
Hinv = [Rt, -Rt*p; 0,0,0,1];
end

