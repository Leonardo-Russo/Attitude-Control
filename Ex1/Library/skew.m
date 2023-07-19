function [A] = skew(v)

% Description: this function creates the skew matrix an input column vector.

n = length(v);

if n ~= 3
    error('Error! The vector provided to the skew() function is not 3x1.\n')
end

h = [];

for i = 1 : n
    h = [h, v];
end

A = cross(h, eye(n));

end