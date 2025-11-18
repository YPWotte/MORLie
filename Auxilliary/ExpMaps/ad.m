function [small_ad] = ad(T)
%AD Small adjoint for twist (w, v) in R^6
w = T(1:3);
v = T(4:6);
sw = skew(w);
small_ad = [sw , zeros(3);...
            skew(v), sw];
end

