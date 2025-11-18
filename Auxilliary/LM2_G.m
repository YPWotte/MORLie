function [G] = LM2_G(Phi,ng,skewg,expG,G0,X,Xnext,Steps,Tol)
%LM_SE3 Levenberg-Marquardt on matrix Lie group G
%   Phi:SE(3) x R^k -> R^k, G0 in G, X in R^k, Xnext in R^k, Steps, Tol
%   Assumes identity inertia to associate algebra gradient with vector
%   Assumes quadratic cost |Phi(H,X) - Xnext|.^2
%
% See also https://en.wikipedia.org/wiki/Levenberg-Marquardt_algorithm

q0 = zeros(ng,1);
dq = 1e-6; lambda = 0;
I = eye(size(G0));

dPhiX = @(G) numJ(@(q)Phi(G*(I+skewg(q)),X),q0, dq);
f = @(G) Phi(G,X) - Xnext;

G = G0;
error = norm(f(G0));
N = 0;
while (error > Tol) && (N < Steps) 
    J = dPhiX(G);
    JJ = J.'*J;
    delta = pinv(JJ+lambda*diag(JJ))*J.'*f(G);
    H_plus = G*expG(-(delta));
    error_plus = norm(f(H_plus));
    if error_plus > error
        lambda = 1 - lambda;
        delta = pinv(JJ+lambda*diag(JJ))*J.'*f(G);
        H_plus = G*expG(-(delta));
        error_plus = norm(f(H_plus));
    end
    if error_plus < error
        N = N+1;
        G = H_plus;
        error = error_plus;
    else
        N = Steps;
    end
end

end

