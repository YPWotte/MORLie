function [H] = LM_SE3(f,H0,Steps,Tol)
%LM_SE3 Levenberg-Marquardt on SE(3)
%   f:SE(3) -> R, H0 in SE(3), Steps, Tol
%   Assumes identity inertia to associate algebra gradient with vector
%
%
% See also https://en.wikipedia.org/wiki/Levenberg-Marquardt_algorithm

q0 = zeros(6,1);
dq = 1e-6; lambda = 0;
df = @(H) numJ(@(q)f(H*(eye4+skew(q))),q0, dq);


H = H0;
error = norm(f(H0));
N = 0;
while (error > Tol) && (N < Steps) 
    J = df(H);
    JJ = J.'*J;
    delta = pinv(JJ+lambda*diag(JJ))*J.'*f(H);
    H_plus = H*exp_se3(-skew(delta));
    error_plus = norm(f(H_plus));
    if error_plus > error
        lambda = 1 - lambda;
        delta = pinv(JJ+lambda*diag(JJ))*J.'*f(H);
        H_plus = H*exp_se3(-skew(delta));
        error_plus = norm(f(H_plus));
    end
    if error_plus < error
        N = N+1;
        H = H_plus;
    else
        N = Steps;
    end
end

end

