function [H] = LM2_SE3(Phi,H0,X,Xnext,Steps,Tol)
%LM_SE3 Levenberg-Marquardt on SE(3)
%   Phi:SE(3) x R^3n -> R^3n, H0 in SE(3), X in R^3n, Xnext in R^3n, Steps, Tol
%   Assumes identity inertia to associate algebra gradient with vector
%   Assumes quadratic cost |Phi(H,X) - Xnext|.^2
%
% See also https://en.wikipedia.org/wiki/Levenberg-Marquardt_algorithm

q0 = zeros(6,1);
dq = 1e-6; lambda = 0;

dPhiX = @(H) numJ(@(q)Phi(H*(eye(4)+skew(q)),X),q0, dq);
f = @(H) Phi(H,X) - Xnext;

H = H0;
error = norm(f(H0));
N = 0;
while (error > Tol) && (N < Steps) 
    J = dPhiX(H);
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

