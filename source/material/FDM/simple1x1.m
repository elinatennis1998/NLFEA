function [ v, dv ] = simple1x1( a, props )
%SIMPLE1X1 Simple material model for scalar/scalar case
%   v = sign(a) * props(1)

v = props(1) * a;
dv = props(1);

%v = props(1) * sign(a);
%dv = 0.0;

end

