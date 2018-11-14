function [tau,intb,vol] = TauE1_3dW(xl,D,nel,lintt,lintl)
% 02/29/2012
% computes tau for edge bubble over a triangle for elasticity problem

tau = zeros(3,3);
ib = 0;
der = 0;
bf = 0;
intb = 0;
vol = 0;
    
lint = lintt*lintl;

for l = 1:lint

    [Wgt,ss] = intpntw(l,lintt,lintl,ib);
    [shl,shld,shls,be] = shlw(ss,6,6,der,bf);
    [shg, shgs, Jdet, be, xs] = shgw(xl,6,shld,shls,6,bf,der,be);

    [b,db] = facebubbleW(ss);
%     db = (db'*inv(xs))';
    bee = (db'/xs)';

    Bmat = [bee(1)  0       0
            0       bee(2)  0
            0       0       bee(3)
            bee(2)  bee(1)  0
            0       bee(3)  bee(2)
            bee(3)  0       bee(1)];

    tau = tau + Jdet*Wgt*Bmat'*D*Bmat;
    intb = intb + Jdet*Wgt*b;
    vol = vol + Jdet*Wgt;

end

tau = inv(tau);
