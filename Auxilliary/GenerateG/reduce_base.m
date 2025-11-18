function [H] = reduce_base(x)
%REDUCE_BASE For vectors in columns x, returns base H of plane orthogonal to x
%   Reduces identity matrix to contain only vectors without components into
%   direction of x, then it uses an SVD on the resulting matrix to return
%   an orthogonal basis.
%   x of zero length are ignored.
n = size(x,1);                                                              % Dimension of space
m = size(x,2);
H = eye(n); % full basis
M = eye(n); % metric
dS = 1e-1;                                                                 % Threshold for singular directions

for i = 1:m
    normx2 = x(:,i)'*M*x(:,i);
    if normx2 > 0
        H = H - (1/normx2)*x(:,i)*(H'*M*x(:,i))';
    end
end

[H,S] = svd(H); % Should use M to return basis orthonormal w.r.t. M
S = sum(diag(S)>dS);                                                         
H = H(:,1:S);

end

