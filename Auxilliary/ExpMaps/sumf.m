function [sum] = sumf(f,th, n)
%SUMF Returns n term series out = sum_i f(i,th) with i from 0 to n
%   f should be a function handle, th should have a defined th^i, n an int 

sum = 0;
for i = 0:n
    sum = sum + f(i,th);
end

end

