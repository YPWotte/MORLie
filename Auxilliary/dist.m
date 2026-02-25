function [d] = dist(S1,S2,pnorm)
%DIST distance between collection of point clouds S1, S2 
% 
% 
%  - S1,S2 are (N*3)*M dimensional arrays, with columns the point cloud at a 
%    given instance, and N the number of points per cloud
%  - d is a 1 by N dimensional array
%
% The distance is the averaged particle distance in the p-norm on R^3, computed as
% d(1,i) = 1/N sum_{n=1}^{N} |S1((n-1)*3+1:n*3,i) - S2((n-1)*3+1:n*3,i)| 

[N3,M] = size(S1);
N = N3/3;
d_Full = zeros(N,M);
for n = 1:N
    d_Full(n,:) = vecnorm(S1((n-1)*3+1:n*3,:) - S2((n-1)*3+1:n*3,:),pnorm); 
end
d = sum(d_Full,1)/N;