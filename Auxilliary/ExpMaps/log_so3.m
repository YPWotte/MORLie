function [wtilde] = log_so3(R)
%LOG_SO3 for input R in SO(3) returns finite wtilde in so(3)
%   Uses specific routine to cover all of SO(3) in a non-closed, non-open set
%   If only R with Tr(R) > -1 are given, then an open cover results.

%% Trace Checking
TrR = trace(R);
A = (R - R')/2; % Antisymmetric part of R
S = (R + R')/2; % Symmetric part of R

%% Main Function
if TrR > -1 && TrR < 3
    %% Generic Case:
    wtilde = acos((TrR-1)/2)/sqrt(-trace(A^2)/2)*A;                         % trace(A^2) could be done cheaper, in principle only 9 multiplications, 8 additions
elseif TrR == 3
    %% Exception if R = eye(3)
    wtilde = A;                                                             % In this case, A = 0
else % i.e. if TrR == -1                                                    % These cases are already double covered by this given chart
    %% Exception if |w| = pi
    w1w2 = sign(S(1,2));
    w1w3 = sign(S(1,3));
    w2w3 = sign(S(2,3));
    w_hat_squared = (S-eye(3))/2;
    if S(1,1)+1 ~= 0  % i.e. w1 =/= 0                                       % Arbitrarily determine which vector to choose for these orientations
        w1 = +sqrt(w_hat_squared(1,1)+1);
        w2 = w1w2*sqrt(w_hat_squared(2,2)+1);
        w3 = w1w3*sqrt(w_hat_squared(3,3)+1);
    elseif S(2,2)+1 ~= 0 % i.e. w2 =/= 0
        w1 = w1w2*sqrt(w_hat_squared(1,1)+1);                               % 0 anyways
        w2 = +sqrt(w_hat_squared(2,2)+1);
        w3 = w2w3*sqrt(w_hat_squared(3,3)+1);
    else % i.e. S(3,3)+1 ~= 0, w3 =/= 0 
        w1 = w1w3*sqrt(w_hat_squared(1,1)+1);
        w2 = w2w3*sqrt(w_hat_squared(2,2)+1);
        w3 = +sqrt(w_hat_squared(3,3)+1);
    end
    wtilde = skew(pi*[w1;w2;w3]);
end

end

