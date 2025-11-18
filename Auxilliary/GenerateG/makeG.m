function [expA,dexpA,dexpA_inv] = makeG(Ei,name,folder)
%MAKEG Generates exponential map, K and Kinv for Lie algebra with basis Ei
%   Functions are saved in ./Generated_name as exp_name, K_name, Kinv_name


if ~exist('folder','var')
    folder = ".";
end

% make folder
dir = sprintf("%s/Generated_%s",folder,name);
mkdir(dir);
adname = sprintf("ad_%s",name);
skewname = sprintf("skew_%s",name);
unskewname = sprintf("unskew_%s",name);
skewdualname = sprintf("skewd_%s",name);
unskewdualname = sprintf("unskewd_%s",name);
expname = sprintf("exp_%s",name);
kname = sprintf("K_%s",name);
kinvname = sprintf("Kinv_%s",name);

% start computing
n = size(Ei,1);
[m,~] = size(Ei{1});
dEi = dual(Ei);
C = structure_coeffs(Ei);

% generatre ad, skew, unskew, dual skew, dual unskew
A = sym("A",[n 1],"real"); %assume(A,"real");
Askew = sym("As",[m m],"real");

ad_fun = simplify(adgl(A,C));
skew_fun = simplify(skewgl(A,Ei));
unskew_fun = simplify(unskewgl(Askew,dEi));
skew_dual_fun = simplify(skewgl(A,dEi));
unskew_dual_fun = simplify(unskewgl(Askew,Ei));

matlabFunction(ad_fun,"File",dir+'/'+adname,"Vars", {A});
matlabFunction(skew_fun,"File",dir+'/'+skewname,"Vars", {A});
matlabFunction(unskew_fun,"File",dir+'/'+unskewname,"Vars", {Askew});
matlabFunction(skew_dual_fun,"File",dir+'/'+skewdualname,"Vars", {A});
matlabFunction(unskew_dual_fun,"File",dir+'/'+unskewdualname,"Vars", {Askew});

% generate exponential map
f = @(x,k) exp(x);
p = simplify(Cayley_Hamilton_Coefficients(f,skewgl(A,Ei)));
pfun = matlabFunction(p,"Vars",{A});
expA_raw = simplify(Cayley_Hamilton_Function(pfun,skewgl(A,Ei),A));
expA = matlabFunction(expA_raw,"File",dir+'/'+expname,"Vars", {A});

% generate derivative of exponential map
x = sym("x");
dexp = (1-exp(-x))/x;
fdexp = @ (x,k) evalinvolved(dexp,x,k);

p = simplify(Cayley_Hamilton_Coefficients(fdexp,adgl(A,C)));
pfun_raw = matlabFunction(p,"Vars",{A});
dexpA_raw = simplify(Cayley_Hamilton_Function(pfun_raw,adgl(A,C),A));
dexpA = matlabFunction(dexpA_raw,"File",dir+'/'+kname,"Vars", {A});


% generate inverse of derivative of exponential map
dexpinv = x/(1-exp(-x));
fdexp_inv = @(x,k) evalinvolved(dexpinv,x,k);
pinv = simplify(Cayley_Hamilton_Coefficients(fdexp_inv,adgl(A,C)));
pfun_rawinv = matlabFunction(pinv,"Vars",{A});
dexpA_inv_raw = simplify(Cayley_Hamilton_Function(pfun_rawinv,adgl(A,C),A));
dexpA_inv = matlabFunction(dexpA_inv_raw,"File",dir+'/'+kinvname,"Vars", {A});

function [out] = evalinvolved(dexp,x,k)
if x == 0
    out = limit(diff(dexp,k),x); 
else
    out = subs(diff(dexp,k),x);
end
end 

end

