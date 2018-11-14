function [tau,intb,Mp,Kp,vol] = Tau9_2d(xl,D,rho,beta,tstep,nel)
% 02/29/2012
% computes tau for edge bubble over a triangle for elasticity problem

tau = zeros(2,2);
Mp = tau;
Kp = tau;
ib = 0;
der = 0;
bf = 1;
intb = 0;
vol = 0;
lint = 9;4;
    
for l = 1:lint

    if nel == 3 || nel == 6
        [Wgt,litr,lits] = intpntt(l,lint,ib);
        [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
        [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
    elseif nel == 4 || nel == 9
        [Wgt,litr,lits] = intpntq(l,lint,ib);
        [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
        [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
    end
    
    db = be(1:2);
    b = be(3);
    Nmat = b*eye(2);
    Bmat = [db(1) 0
            0 db(2)
            db(2) db(1)];

    Mp = Mp + Jdet*Wgt*(Nmat'*rho*Nmat);
    Kp = Kp + Jdet*Wgt*(Bmat'*D*Bmat);
    intb = intb + Jdet*Wgt*b;
    vol = vol + Jdet*Wgt;

end

tau = tau + (1/(beta*tstep^2)*Mp+Kp);

tau = inv(tau);
