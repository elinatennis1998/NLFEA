function [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,a,b,c,rho,beta,tstep,nel)
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

Bumat = zeros(3,nel*2);
Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];
    
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
    
            
    % Form B matrix
    for ie = 1:nel

      Bumat(Bcol1,(ie-1)*2+1) = shg(ie,col1);
      Bumat(Bcol2,(ie-1)*2+2) = shg(ie,col2);

    end
            
    E = (Bumat*ulres);
    E(3) = E(3)/2; %for extra factor of 2 on shear strain
    Ekk = E(1) + E(2);
    
    db = be(1:2);
    bub = be(3);
    Nmat = bub*eye(2);
    Bmat = [db(1) 0
            0 db(2)
            db(2) db(1)];
    cmat = [a+b+2*c*E(1) a c*E(3) 
            a a+b+2*c*E(2) c*E(3)
            c*E(3) c*E(3)  b/2+c/2*Ekk];

    Mp = Mp + Jdet*Wgt*(Nmat'*rho*Nmat);
    Kp = Kp + Jdet*Wgt*(Bmat'*cmat*Bmat);
    intb = intb + Jdet*Wgt*bub;
    vol = vol + Jdet*Wgt;

end

tau = tau + (1/(beta*tstep^2)*Mp+Kp);

tau = inv(tau);
%tau = zeros(2,2);
