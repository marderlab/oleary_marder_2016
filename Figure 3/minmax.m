function [idx, m] = minmax(x, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = length(y);

if length(x) ~= N
    error('x and y must be equal length');
end

for i=1:N
    if find(x(i) == y,1,'first') <= i
        idx = i;
        m = x(i);
        return;
    end
end

end
