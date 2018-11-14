function [b,bd] = bubbleCST(c1,c2,xs)

bd = zeros(1,2);
twent7 = 27.d0;
c3 = 1.d0 - c1 - c2;

bd(1) = twent7*c2*(c3-c1);
bd(2) = twent7*c1*(c3-c2);
b = twent7*c1*c2*c3;

bd = bd*xs;