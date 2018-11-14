function [dvec,dddt,dmax,dstate,bonddbond2,dbond] = EllipseCZM(tvtr,jumpu,dmax,dvec,nvec,sigmax,beta,Hc,dc,rp,dstate,elast_initial,dtol,ndm)
%
% 08/10/2014
% Tim Truster
%
% Debonding interface model, unequal tension and shear limit, non-penetration
% contact. Implementation borrowed from NL_Elem4_3d in Interface4 folder.
%
% elast_initial = 1 for the first iteration of a new load step, for example
% determined by transient == 0 and initia == 1; studies have shown that
% initializing the stiffness matrix with the damage tangent yields stable
% results compared to using the elastic/unloading tangent.
% dstate is a flag for what portion of the CZM model is active, e.g. damage
% or adhesion.

beta2 = beta^2;
            
sn = nvec'*tvtr;                              %normal stress
un = jumpu'*nvec;                             %normal gap
tn = sn+rp*un;                                %normal traction
tss = tvtr + rp*jumpu - tn*nvec;              %shear traction v
tt = sqrt(tss'*tss);                          %shear traction

if tn >= 0 % tension

    tvec = tvtr + rp*jumpu;
    tb = sqrt(tn^2+beta2*tt^2);

else % compression

    tvec = tss;
    tb = beta*tt;

end

normtvec = sqrt(tvec'*tvec);
if dmax >= dc
    psik = 0;
else
    psik = sigmax - Hc*dmax;
end

%%            
if tn >= 0 % tension

  if dmax == 0 && sqrt(tn^2+tt^2/beta2) <= sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dmax = 0;
    bonddbond2 = 0;
    dstate = 1;

  elseif dmax == 0 && sqrt(tn^2+tt^2/beta2) > sigmax % damage from initial adhesion

    [dvec,dtdd,deq] = damageNR(tvec,nvec,dvec,normtvec,sigmax,dc,rp,beta2,ndm);
    if deq >= dc % opening/sliding
        dvec = tvec/rp;
        dddt = 1/rp*eye(ndm);
        dn = dvec'*nvec;
        dss = dvec - dn*nvec;
        dt = sqrt(dss'*dss);
        deq = sqrt(dn^2+beta2*dt^2);
    else
        dddt = inv(dtdd);
    end
    dmax = max(dmax,deq);
    bonddbond2 = -1;
    dstate = 4;

  elseif dmax > 0 && tb < rp*dmax % crumpling

    dvec = tvec/rp;
    dddt = 1/rp*eye(ndm);
    dn = dvec'*nvec;
    dss = dvec - dn*nvec;
    dt = sqrt(dss'*dss);
    deq = sqrt(dn^2+beta2*dt^2);
    dmax = max(dmax,deq);
    bonddbond2 = -1;
    dstate = 2;

  elseif dmax > 0 && rp*dmax <= normtvec && ... %unloading
          (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 <= 1 && ~(dstate == 4 && elast_initial)) || ...
          (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 <= 1-dtol && (dstate == 4 && elast_initial)) % make it harder to unload during initialization step
  
    [dvec,dddt,deq] = unloadingNR(tn,tt,tss,nvec,dmax,sigmax,dc,rp,beta2,ndm);
    bonddbond2 = -1;
    dstate = 3;
    
  elseif ... % damage
      (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 > 1 && ~(dstate == 4 && elast_initial)) || ...
      (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 > 1-dtol && (dstate == 4 && elast_initial)) % damage predictor

    [dvec,dtdd,deq] = damageNR(tvec,nvec,dvec,normtvec,sigmax,dc,rp,beta2,ndm);
    if deq >= dc % opening/sliding
        dvec = tvec/rp;
        dddt = 1/rp*eye(ndm);
        dn = dvec'*nvec;
        dss = dvec - dn*nvec;
        dt = sqrt(dss'*dss);
        deq = sqrt(dn^2+beta2*dt^2);
    else
        dddt = inv(dtdd);
    end
    dmax = max(dmax,deq);
    bonddbond2 = -1;
    dstate = 4;

  end

else % compression

  if dmax == 0 && tt/beta <= sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dmax = 0;
    bonddbond2 = 0;
    dstate = 1;

  elseif dmax > 0 && tb < rp*dmax % crumpling

    dvec = tvec/rp;
    dddt = 1/rp*(eye(ndm) - nvec*nvec');
    deq = beta*sqrt(dvec'*dvec);
    dmax = max(dmax,deq);
    bonddbond2 = -1;
    dstate = 2;
  elseif dmax > 0 && rp*dmax <= tb && ... %unloading
          (tb <= rp*dmax + beta2*psik && ~(dstate == 4 && elast_initial)) || ...
          (tb <= rp*dmax + beta2*psik - dtol && (dstate == 4 && elast_initial)) % make it harder to unload during initialization step

    ntilda = tss/tt;
    dvec = (dmax/beta)*ntilda;
    dddt = (dmax/beta)/tt*(eye(ndm) - ntilda*ntilda' - nvec*nvec');
    bonddbond2 = -1;
    dstate = 3;

  elseif ... % damage
      (tb > rp*dmax + beta2*psik && ~(dstate == 4 && elast_initial)) || ...
      (tb > rp*dmax + beta2*psik - dtol && (dstate == 4 && elast_initial)) % damage predictor

    deq = beta*(tt - beta*sigmax)/(rp - beta2*Hc); % assume only damage
    ntilda = tss/tt;
    if deq >= dc % opening/sliding
        dvec = tvec/rp;
        deq = beta*sqrt(dvec'*dvec);
        dddt = 1/rp*(eye(ndm) - nvec*nvec');
    else % only damage
        dvec = (deq/beta)*ntilda;
        dddt = 1/(1*(rp - beta2*Hc))*(ntilda*ntilda') + (deq/beta)/tt*(eye(ndm) - ntilda*ntilda' - nvec*nvec');
    end
    dmax = deq;
    bonddbond2 = -1;
    dstate = 4;

  end

end

if dmax > 0 && tn >= 0
    if dmax > dc
        dbond = 2;
    else
        dbond = 1;
    end
elseif dmax > 0 && tn < 0
    if dmax > dc
        dbond = -2;
    else
        dbond = -1;
    end
else
    dbond = 0;
end