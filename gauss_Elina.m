function [r,w] = gauss_Elina(lint)
%Elina Geut. Last modified 11/7/2018
%   Function that conatins weight functions ang integration points

w = zeros(1,5);
if lint == 1
    r(1) = -1;
    r(2) = 1;
    w(1) = 2;

elseif lint == 2
    r(1) = -1/(sqrt(3));
    r(2) = 1/(sqrt(3));
    w(1) = 1;
    w(2) = 1;

elseif lint == 3
    r(1) = -0.7745966692;
    r(2) = 0;
    r(3) = 0.7745966692;
    w(1) = 5/9;
    w(2) = 8/9;
    w(3) = 5/9;
elseif lint == 4
    r(1) = -0.8611363116;
    r(2) = -0.3399810436;
    r(3) = 0.3399810436;
    r(4) = 0.8611363116;
    w(1) = 0.3478548451;
    w(2) = 0.6521451549;
    w(3) = w(2);
    w(4) = w(1);
elseif lint == 5
    r(1) = -0.9061798459;
    r(2) = -0.5384693101;
    r(3) = 0;
    r(4) = 0.5384693101;
    r(5) = 0.9061798459;
    w(1) = 0.2369268851;
    w(2) = 0.4786286705;
    w(3) = 0.5688888889;
    w(4) = w(2);
    w(5) = w(1);
end
end

