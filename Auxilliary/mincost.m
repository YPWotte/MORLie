function [cost_value] = mincost(A,X,X_A,F,x)
%MINCOST Summary of this function goes here
%   Detailed explanation goes here
cost = @(B) F(x,X_A(B,x)-X(x));

smin = fminsearch(@(s)cost(s*A),0);
cost_value = cost(smin*A);
end

