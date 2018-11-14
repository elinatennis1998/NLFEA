function [dvec,dddt,deq] = unloadingNR(tn,tt,tss,nvec,k,sc,dc,rp,beta2,ndm)
%
% 06/26/2012
% Tim Truster
%
% Function to compute updated delta vector by Newton-Raphson for unloading
%
% Updated 07/13/2012 to include quartic equation solver and a check for
% initial guess of rho = psik/k

if k >= dc

    dvec = (tn*nvec + tss)/rp;
    dddt = 1/rp*eye(ndm);
    dn = dvec'*nvec;
    dss = dvec - dn*nvec;
    dt = sqrt(dss'*dss);
    deq = sqrt(dn^2+beta2*dt^2);

else

residratio = 1e-11;
itermax = 10;

I3 = eye(ndm);
psik = sc*(1 - k/dc);

r2 = psik/k;
f2 = k^2 - tn^2/(rp + r2)^2 - beta2*tt^2/(rp + beta2*r2)^2;

if abs(f2) < residratio
    
    rho = r2;
    KNR = -2*tn^2/(rp + rho)^3 - 2*beta2^2*tt^2/(rp + beta2*rho)^3;
    exitflag = 1;
    
else

% Initialize
% rho = 0.95*psik/k;
options = optimset('TolFun',1e-12,'Display','off','Jacobian','on');%
r1 = 0;
f1 = k^2 - tn^2/(rp + r1)^2 - beta2*tt^2/(rp + beta2*r1)^2;
r2 = psik/k;
f2 = k^2 - tn^2/(rp + r2)^2 - beta2*tt^2/(rp + beta2*r2)^2;
m = (f2-f1)/(r2-r1);
rho0 = r2 - f2/m;
[rho,fval,exitflag] = fsolve(@(rho) unl(k,rp,beta2,tn,tt,rho),rho0,options);
KNR = -2*tn^2/(rp + rho)^3 - 2*beta2^2*tt^2/(rp + beta2*rho)^3;

end

dvec = tn/(rp + rho)*nvec + tss/(rp + beta2*rho);
dn = dvec'*nvec;
dss = dvec - dn*nvec;
dt = sqrt(dss'*dss);
deq = sqrt(dn^2+beta2*dt^2);

nn = (nvec*nvec');
mvec = tss/tt;
term1 = 1/(rp + rho)*nn;
term2 = 1/(rp + beta2*rho)*(I3 - nn);
numern = tn/(rp + rho)^2;
numert = beta2*tt/(rp + beta2*rho)^2;
term3 = (-2*numern*numern/KNR)*nn + numern*(-2*numert/KNR)*(nvec*mvec');
term4 = (-2*numert*numern/KNR)*mvec*nvec' + (-2*numert*numert/KNR)*(mvec*mvec');
dddt = term1 + term2 - term3 - term4;

if exitflag < 1
    disp('unloadingNR')
    exitflag
end

end