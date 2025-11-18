addpath(genpath('../ExpMaps'));

name = "SE3";
Ei = cell(6,1);
Ei_flat = eye(6);
for i = 1:6
    Ei{i} = skew(Ei_flat(:,i));
end
makeG(Ei,name);

% dEi = dual(Ei);
% 
% % dEi1 = reshape(dEi{1},[],1);
% % Ei1 = reshape(Ei{1},[],1);
% % dEi1.'*Ei1;
% 
% T = rand(6,1);
% Tskew = skewgl(T,Ei);
% zero = T - unskewgl(Tskew,dEi);
% 
% C = structure_coeffs(Ei);
% 
% % generate exponential map
% syms A [6 1] real;
% f = @(x,k) exp(x);
% p = Cayley_Hamilton_Coefficients(f,skewgl(A,Ei));
% pfun = matlabFunction(real(p),"Vars",{A});
% expA = @(skewA) Cayley_Hamilton_Function(pfun,skewA,unskewgl(skewA,dEi));
% 
% % generate derivative of exponential map
% syms x;
% dexp = (1-exp(-x))/x;
% fdexp = @(x,k) limit(diff(dexp,k),x);
% p = simplify(Cayley_Hamilton_Coefficients(fdexp,adgl(A,C)));
% pfun_raw = matlabFunction(real(p),"Vars",{A});
% dexpA = @(A) Cayley_Hamilton_Function(pfun_raw,adgl(A,C),A);
% 
% 
% % generate inverse of derivative of exponential map
% dexpinv = x/(1-exp(-x));
% fdexp = @(x,k) limit(diff(dexpinv,k),x);
% pinv = simplify(Cayley_Hamilton_Coefficients(fdexp,adgl(A,C)));
% pfun_rawinv = matlabFunction(real(pinv),"Vars",{A});
% dexpA_inv = @(A) Cayley_Hamilton_Function(pfun_rawinv,adgl(A,C),A);



% simplify(eig(fA(skewgl(A,Ei)))
% syms x;
% logx = log(x);
% flog = @(x,k) double(subs(diff(logx,k),x));
% 
% p = Cayley_Hamilton_Coefficients(f,skewgl(A,Ei));
% pfun_raw = matlabFunction(real(p),"Vars",{A});
% pfun = @(skewA) pfun_raw(unskewgl(skewA,dEi));
% fA = @(skewA) Cayley_Hamilton_Function(pfun,skewA);



% zero2 = expA(Tskew) - expm(Tskew);
%pfun = matlabFunction(p,"File","matrixF","Vars",{T});

