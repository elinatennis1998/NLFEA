function [dvec,KNR,deq] = damageNR(tvec,nvec,dvecn_1,normtvec,sc,dc,rp,beta2,ndm)
%
% 06/26/2012
% Tim Truster
%
% Function to compute updated delta vector by Newton-Raphson for damage

residratio = 1e-12;
itermax = 10;

I3 = eye(ndm);

% Initialize
dvec = dvecn_1;
if norm(dvec)/dc < residratio
    tn = tvec'*nvec;
    tss = tvec - tn*nvec;
    tt = sqrt(tss'*tss);
    dvec = 1e-7*dc*(tn*nvec + tss/beta2)/sqrt(tn^2+tt^2/beta2);
end
dn = dvec'*nvec;
dss = dvec - dn*nvec;
dt = sqrt(dss'*dss);
deq = sqrt(dn^2+beta2*dt^2);
dhat = dn*nvec + beta2*dss;
% if deq > 0
resid = tvec - rp*dvec - sc*(1-deq/dc)*(dhat/deq);
KNR = rp*I3 + sc*(1/deq - 1/dc)*(nvec*nvec' + beta2*(I3 - nvec*nvec')) - (sc/deq)*(dhat/deq)*(dhat/deq)';
% else
% % use beta = 1 for simplicity
% ntilda = tvec/normtvec;
% resid = tvec - rp*dvec - sc*ntilda;
% KNR = rp*I3 + (sc*(1-deq/dc)
% end
rnorm = norm(resid);
if rnorm < 50*residratio %already converged, i.e. linear
    residtol = 1.2*rnorm;
else
    residtol = max(residratio*rnorm,40*eps);
end
iter = 0;

while iter < itermax && rnorm > residtol
   
    iter = iter + 1;
    
    ddvec = KNR\resid;
    
    dvec = dvec + ddvec;
    if norm(dvec)/dc < residratio
        return
    end
    dn = dvec'*nvec;
    dss = dvec - dn*nvec;
    dt = sqrt(dss'*dss);
    deq = sqrt(dn^2+beta2*dt^2);
    dhat = dn*nvec + beta2*dss;
    resid = tvec - rp*dvec - sc*(1-deq/dc)*(dhat/deq);
    KNR = rp*I3 + sc*(1/deq - 1/dc)*(nvec*nvec' + beta2*(I3 - nvec*nvec')) - (sc/deq)*(dhat/deq)*(dhat/deq)';
    rnorm = norm(resid);
    
end
if iter > 5
    disp('damageNR')
    iter
end