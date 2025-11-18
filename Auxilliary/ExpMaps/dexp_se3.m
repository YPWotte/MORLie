function [K] = dexp_se3(q,type)
%DEXP_SE3 Calculates derivative of exponential map using alternative method
%   Maths found in "Derivative of the Expoential Map, Ethan Eade 2018"
%   type defaults to 'normal', but 'inverse' may also be given to get K^(-1)
%% Check additional inputs
if ~exist('type','var')
    type = 'normal';
end

%% Extract w, v and angle theta
w = q(1:3);
v = q(4:6);

th = sqrt(w'*w);

%% Calculation of coefficients
if th > 0.01
    a = sin(th)/th;
    b = (1-cos(th))/th^2;
    c = (1 - a)/th^2;
    f = (a-2*b)/th^2;
    g = (b-3*c)/th^2;
    if (abs(sin(th)) < 0.1) % equivalent formulas, but different singularities
    e = (b - a/2)/(1-cos(th));
    else
        e = (b-2*c)/(2*a);
    end
else    % Series expansions to deal with ill numerical conditioning
    n = max(4,-ceil(20*log(10)/log(th)));                                   % keeps relative error around th^5, absolute error around 10^-20
    as = @(i,th) (-1)^i/factorial(2*i+1)*th^(2*i);
    bs = @(i,th) (-1)^i/factorial(2*i+2)*th^(2*i);
    cs = @(i,th) (-1)^i/factorial(2*i+3)*th^(2*i);
    fs = @(i,th) ((-1)^(i+1)/factorial(2*i+3)-2*(-1)^(i+1)/factorial(2*i+4))*th^(2*i);
    gs = @(i,th) ((-1)^(i+1)/factorial(2*i+4)-3*(-1)^(i+1)/factorial(2*i+5))*th^(2*i);
    a = sumf(as,th,n); % 1;
    b = sumf(bs,th,n); % 1/2;
    c = sumf(cs,th,n); % 1/6;
    f = sumf(fs,th,n); % -1/12;
    g = sumf(gs,th,n); % -1/60;
    e = -f/(2*b);      % 1/12;
end

Q = f*skew(w) + g*skew(w)^2;
W = Q - 2*c*eye(3);

%% Calculate output
if type == "inverse"
    B = b*skew(v)+c*(w*v'+v*w')+(w'*v)*W;
    Dw_inv = eye(3) - skew(w)/2 + e *skew(w)^2;
    K = [Dw_inv, -Dw_inv*B*Dw_inv; zeros(3,3), Dw_inv]'; 
else
    Dw = a*eye(3) + b*skew(w) + c*w*w';
    Dv = b*skew(v) + c*(w*v'+v*w')+(w'*v)*W;
    K = [Dw, Dv; zeros(3,3), Dw]';
end

end