function [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter,tvec] = RadialCZME(tvtr,jumpu,dmax,dvec, ...
    nvec,sigmax,Hc,dc,rp,dstate,iter,elast_initial,itchat,chatter,bonddbond1,dtol,ndm)
%
% 12/20/2015
% Tim Truster
%
% Debonding interface model, equal tension and shear limit, non-penetration
% contact. Termed radial because it is like radial return.
%
% elast_initial = 1 for the first iteration of a new load step, for example
% determined by transient == 0 and initia == 1; studies have shown that
% initializing the stiffness matrix with the damage tangent yields stable
% results compared to using the elastic/unloading tangent.
% dstate is a flag for what portion of the CZM model is active, e.g. damage
% or adhesion.
%
% NOTE: currently, only the main portions of the material model have been
% updated, not the chatter portions.
%
% Added elastic slope for unloading

            
normsig = nvec'*tvtr;
dn = jumpu'*nvec;
tn = normsig+rp*dn;

if tn >= 0 % tension

    tvec = tvtr + rp*jumpu;

else % compression

    tvec = tvtr + rp*jumpu - tn*nvec;

end

normtvec = sqrt(tvec'*tvec);
if dmax >= dc
    psik = 0;
else
    psik = sigmax - Hc*dmax;
end  

%% Branch on region of potential function
if tn >= 0 % tension

 if iter > itchat && chatter == 1

  if dmax == 0 && normtvec <= 1.5*sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dddu = rp*dddt;
    dmax = 0;
    dstate = 1;
% 
%   elseif dmax > 0 && normtvec < rp*dmax % crumpling
% 
%     dvec = tvec/rp;
%     dddt = 1/rp*eye(ndm);
%     dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
%     dstate = 2;

  elseif dmax > 0 && normtvec < rp*dmax + psik && (dstate ~= 4 && iter == 0) %unloading

    Kc = (sigmax-sigmax*dmax/dc)/dmax;
    dvec = 1/(Kc+rp)*tvec;
    dddt = 1/(Kc+rp)*(eye(ndm));
    dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
    dstate = 3;
                
  elseif rp*dmax + psik <= normtvec || (dstate == 4 && iter == 0) % damage

    ntilda = tvec/normtvec;
    deq = (normtvec - sigmax)/(rp - Hc); % assume only damage
    if deq >= dc % opening/sliding
        deq = normtvec/rp;
        dddt = 1/rp*eye(ndm);
        dddu = rp*dddt;
    else % only damage
        dddt = 1/(rp - Hc)*(ntilda*ntilda') + deq/normtvec*(eye(ndm) - ntilda*ntilda');
        dddu = rp*dddt;
    end
    dvec = deq*ntilda;
    dmax = deq;
    dstate = 4;

  end

  
 else
         
     
  if dmax == 0 && normtvec <= sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dddu = rp*dddt;
    dmax = 0;
    bonddbond2 = 0;
    dstate = 1;
% 
%   elseif dmax > 0 && normtvec < rp*dmax % crumpling
% 
%     dvec = tvec/rp;
%     dddt = 1/rp*eye(ndm);
%     dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
%     bonddbond2 = -1;
%     dstate = 2;

  elseif dmax > 0 && ... %unloading
          (normtvec <= rp*dmax + psik && ~(dstate == 4 && elast_initial)) || ...
          (normtvec <= rp*dmax + psik - dtol && (dstate == 4 && elast_initial)) % make it harder to unload during initialization step

    Kc = (sigmax-sigmax*dmax/dc)/dmax;
    dvec = 1/(Kc+rp)*tvec;
    dddt = 1/(Kc+rp)*(eye(ndm));
    dddu = rp*dddt;
%                 deq = sqrt(dvec'*dvec);
%                 dmax = max(dmax,deq);
    bonddbond2 = -1;
    dstate = 3;

  elseif ... % damage
      (rp*dmax + psik < normtvec && ~(dstate == 4 && elast_initial)) || ...
      (rp*dmax + psik - dtol < normtvec && (dstate == 4 && elast_initial)) % damage predictor

    ntilda = tvec/normtvec;
    deq = (normtvec - sigmax)/(rp - Hc); % assume only damage
    if deq >= dc % opening/sliding
        deq = normtvec/rp;
        dddt = 1/rp*eye(ndm);
        dddu = rp*dddt;
    else % only damage
        dddt = 1/(rp - Hc)*(ntilda*ntilda') + deq/normtvec*(eye(ndm) - ntilda*ntilda');
        dddu = rp*dddt;
    end
    dvec = deq*ntilda;
    dmax = deq;
    bonddbond2 = -1;
    dstate = 4;
    if dmax < 0
        dmax
    end

  end

 end

 
else % compression

    
 if iter > itchat && chatter == 1

  if dmax == 0 && normtvec <= 1.5*sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dddu = rp*dddt;
    dmax = 0;
    dstate = 1;
% 
%   elseif dmax > 0 && normtvec < rp*dmax % crumpling
% 
%     dvec = tvec/rp;
%     dddt = 1/rp*(eye(ndm) - nvec*nvec');
%     dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
%     dstate = 2;

  elseif dmax > 0 && normtvec < rp*dmax + psik && (dstate ~= 4 && iter == 0) %unloading

    Kc = (sigmax-sigmax*dmax/dc)/dmax;
    dvec = 1/(Kc+rp)*tvec;
    dddt = 1/(Kc+rp)*(eye(ndm) - nvec*nvec');
    dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
    if norm(dddt - dddt',1) > 1e-9
        keyboard
    end
    dstate = 3;

  elseif rp*dmax + psik <= normtvec || (dstate == 4 && iter == 0) % damage

    ntilda = tvec/normtvec;
    deq = (normtvec - sigmax)/(rp - Hc); % assume only damage
    if deq >= dc % opening/sliding
        deq = normtvec/rp;
        dddt = 1/rp*(eye(ndm) - nvec*nvec');
        dddu = rp*dddt;
    else % only damage
        dddt = 1/(rp - Hc)*(ntilda*ntilda') + deq/normtvec*(eye(ndm) - ntilda*ntilda' - nvec*nvec');
        dddu = rp*dddt;
    end
    dvec = deq*ntilda;
    dmax = deq;
    if norm(dddt - dddt',1) > 1e-9
        keyboard
    end
    dstate = 4;

  end

  
 else
     
     
  if dmax == 0 && normtvec <= sigmax % perfect adhesion

    dvec = zeros(ndm,1);
    dddt = zeros(ndm,ndm);
    dddu = rp*dddt;
    dmax = 0;
    bonddbond2 = 0;
    dstate = 1;
% 
%   elseif dmax > 0 && normtvec < rp*dmax % crumpling
% 
%     dvec = tvec/rp;
%     dddt = 1/rp*(eye(ndm) - nvec*nvec');
%     dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
%     bonddbond2 = -1;
%     dstate = 2;
    
  elseif dmax > 0 && ... %unloading
          (normtvec <= rp*dmax + psik && ~(dstate == 4 && elast_initial)) || ...
          (normtvec <= rp*dmax + psik - dtol && (dstate == 4 && elast_initial)) % make it harder to unload during initialization step

    Kc = (sigmax-sigmax*dmax/dc)/dmax;
    dvec = 1/(Kc+rp)*tvec;
    dddt = 1/(Kc+rp)*(eye(ndm) - nvec*nvec');
    dddu = rp*dddt;
%     deq = sqrt(dvec'*dvec);
%     dmax = max(dmax,deq);
    if norm(dddt - dddt',1) > 1e-9
        keyboard
    end
    bonddbond2 = -1;
    dstate = 3;

  elseif ... % damage
      (rp*dmax + psik < normtvec && ~(dstate == 4 && elast_initial)) || ...
      (rp*dmax + psik - dtol < normtvec && (dstate == 4 && elast_initial)) % damage predictor

    ntilda = tvec/normtvec;
    deq = (normtvec - sigmax)/(rp - Hc); % assume only damage
    if deq >= dc % opening/sliding
        deq = normtvec/rp;
        dddt = 1/rp*(eye(ndm) - nvec*nvec');
        dddu = rp*dddt;
    else % only damage
        dddt = 1/(rp - Hc)*(ntilda*ntilda') + deq/normtvec*(eye(ndm) - ntilda*ntilda' - nvec*nvec');
        dddu = rp*dddt;
    end
    dvec = deq*ntilda;
    dmax = deq;
    if norm(dddt - dddt',1) > 1e-9
        keyboard
    end
    bonddbond2 = -1;
    dstate = 4;
    if dmax < 0
        dmax
    end

  end

 end

end

tvec(1) = tvec(1) - rp*dvec(1);
tvec(2) = tvec(2) - rp*dvec(2);
if (ndm == 3) % then
tvec(3) = tvec(3) - rp*dvec(3);
end

%% Determine flags for debonding
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
if iter <= itchat
  if sign(bonddbond1) ~= sign(bonddbond2)
    chatter = 1;
  else
    chatter = 0;
  end
end