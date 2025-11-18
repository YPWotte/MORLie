function [h] = log_se3(H)
%LOG_SE3 for input H in SE(3) returns corresponding finite twist
%   Uses specific routine to cover all of SE(3) in a non-closed, non-open set
%   If only H with Tr(H) > 0 are given, then an open cover results.

%% Trace Checking
TrH = trace(H);

%% Main Function
R = H(1:3,1:3);
wtilde = log_so3(R);
normw = norm(unskew(wtilde));
if TrH > 0 && TrH < 4
    Scale = (2*sin(normw)-normw*(1+cos(normw)))/(2*normw^2*sin(normw));
    AInv = eye(3) - wtilde/2 + Scale*wtilde^2;
    v = AInv*H(1:3,4);
    h = [wtilde, v; [0,0,0,0]];
elseif TrH == 4
    v = H(1:3,4);
    h = [wtilde, v; [0,0,0,0]];
else
    Scale = 1/(normw^2);
    AInv = eye(3) - wtilde/2 + Scale*wtilde^2;
    v = AInv*H(1:3,4);
    h = [wtilde, v; [0,0,0,0]];
end

% A = eye(3) + wtilde/normw^2*(1-cos(normw))+ wtilde^2/normw^3*(normw-sin(normw));
end
