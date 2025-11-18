function [BigAd] = AdH(H)
%ADH Returns Adjoint representation of 4x4 input H in SE(3)
%   Assumes that twists are given rotation first (T = (w,v))
R = H(1:3,1:3);
p = H(1:3,4);

BigAd = zeros(6);
BigAd(1:3,1:3) = R;
BigAd(4:6,1:3) = skew(p)*R;
BigAd(4:6,4:6) = R;

end

