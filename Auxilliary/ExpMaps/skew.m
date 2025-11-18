function [Skew_x] = skew(x)
%SKEW Returns skew symmetric matrix for x in R^3 or tilde form for x in R^6
%   
if length(x) == 3
    Skew_x = [ 0, -x(3), x(2); ...
           x(3), 0 , -x(1);...
           -x(2), x(1), 0];
elseif length(x) == 6
    Skew_x = [skew(x(1:3)), x(4:6); zeros(1,4)];
end

