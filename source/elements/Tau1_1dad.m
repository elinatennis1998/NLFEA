function [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,D,rho,coeffm,coeffk,nel,ndm,nen)
% 02/29/2012
% computes tau for edge bubble over 1d linear element

tau = 0;
Mp = tau;
Kp = tau;
ib = 0;
der = 0;
bf = 1;
intb = 0;
vol = 0;
lint = 3;2; % Important to use full integration for tau

sw = int1d(lint);
[shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
    
for l = 1:lint
        
    Wgt = sw(2,l);
    [shg, shgs, Jdet, be(l,:)] = shg1d(xl,ndm,nel,shld(l,:),shls(l,:),nen,bf,der,be(l,:));
    
    db = be(l,2);
    b = be(l,1);
    Nmat = b;
    Bmat = db;

    Mp = Mp + Jdet*Wgt*(Nmat'*rho*Nmat);
    Kp = Kp + Jdet*Wgt*(Bmat'*D*Bmat);
    intb = intb + Jdet*Wgt*b;
    vol = vol + Jdet*Wgt;

end

tau = tau + (coeffm*Mp+coeffk*Kp);

tau = inv(tau);
