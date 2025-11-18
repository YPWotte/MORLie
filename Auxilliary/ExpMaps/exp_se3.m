function [exph] = exp_se3(h)
%EXPH For input h in se(3) returns matrix exponential
%   Uses matrix exponential exp_so3 by Rodrigues' formula for w in so(3)

w_tilde = h(1:3,1:3);
v = h(1:3,4);
w = unskew(w_tilde);
normw = norm(w);

expw = exp_so3(w_tilde);
exph = zeros(4);

exph(1:3,1:3) = expw;

if normw > 0
    exph(1:3,4) = (1/normw^2)*((eye(3)-expw)*skew(w)*v+(w*w')*v);
else
    exph(1:3,4) = v;
end

exph(4,4) = 1;

end

