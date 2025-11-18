addpath(genpath('../ExpMaps'));

n = 3; %m = 3;
name = "SE2";
Ei = cell(n,1);
%Ei_flat = eye(m);
Ei{1} = [0, -1, 0; 1, 0, 0; 0,0,0];
Ei{2} = [0, 0, 1; 0, 0, 0; 0,0,0];
Ei{3} = [0, 0, 0; 0, 0, 1; 0,0,0];
makeG(Ei,name);


% dEi = dual(Ei);

% dEi1 = reshape(dEi{1},[],1);
% Ei1 = reshape(Ei{1},[],1);
% dEi1.'*Ei1;

% T = rand(n,1);
% Tskew = skewgl(T,Ei);
% zero = T - unskewgl(Tskew,dEi);
% 
% C = structure_coeffs(Ei);

% % generate exponential map
% syms A [n 1] real;
% f = @(x,k) exp(x);
% p = simplify(Cayley_Hamilton_Coefficients(f,skewgl(A,Ei)));
% pfun = matlabFunction(real(p),"Vars",{A});
% expA_raw = simplify(Cayley_Hamilton_Function(pfun,skewgl(A,Ei),A));
% expA = matlabFunction(expA_raw,"Vars", {A});
% 
% 
% % generate derivative of exponential map
% syms x real;
% dexp = (1-exp(-x))/x;
% fdexp = @ (x,k) evalinvolved(dexp,x,k);
% 
% p = simplify(Cayley_Hamilton_Coefficients(fdexp,adgl(A,C)));
% pfun_raw = matlabFunction(real(p),"Vars",{A});
% dexpA_raw = simplify(Cayley_Hamilton_Function(pfun_raw,adgl(A,C),A));
% dexpA = matlabFunction(dexpA_raw,"Vars", {A});
% 
% 
% % generate inverse of derivative of exponential map
% dexpinv = x/(1-exp(-x));
% fdexp_inv = @(x,k) limit(diff(dexpinv,k),x);
% pinv = simplify(Cayley_Hamilton_Coefficients(fdexp_inv,adgl(A,C)));
% pfun_rawinv = matlabFunction(real(pinv),"Vars",{A});
% dexpA_inv = @(A) Cayley_Hamilton_Function(pfun_rawinv,adgl(A,C),A);
% 
% fdexp_inv = @(x,k) limit(diff(dexp,k),x); fdexp2 = @(x,k) subs(diff(dexp,k),x);
% 
% 
% function [out] = evalinvolved(dexp,x,k)
% if x == 0
%     out = limit(diff(dexp,k),x); 
% else
%     out = subs(diff(dexp,k),x);
% end
% end