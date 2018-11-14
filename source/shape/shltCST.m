function shl = shltCST(c1,c2)


shl = zeros(3,1);

c3 = 1.0d0 - c1 - c2;

shl(2) = c1;
shl(3) = c2;
shl(1) = c3;
