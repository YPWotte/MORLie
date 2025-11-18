function [exp_w_tilde] = exp_so3(w_tilde)
%EXP_SO3 for w_tilde in se(3), returns matrix exponential using Rodrigues'
%formula

w = unskew(w_tilde);
normtheta = sqrt(w'*w);
if normtheta > 0 
    unit_w_tilde = w_tilde/normtheta;
    exp_w_tilde = eye(3) + unit_w_tilde*sin(normtheta)...
                         + unit_w_tilde^2*(1-cos(normtheta));
else % i.e. if normtheta = 0
    exp_w_tilde = eye(3);
end
end

