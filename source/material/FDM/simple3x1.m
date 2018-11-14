function [ v, dv ] = simple3x1(x, a, props )
%SIMPLE1X3 Simple material model for vector/scalar case
% Matches weird thing that works on my python version
% x: position
% a: field
% props: properties

v0 = props(1);
c = props(2);

if x < c
    v = v0;
elseif x > c
    v = -v0;
else
    v = 0.0;
end

dv = zeros(3,1);

end

