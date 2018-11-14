function [zeta,dzdt,kappa,dstate,bonddbond2,dbond,chatter,tvec] = GeneralCZM(tvtr,jumpu,kappa,zeta, ...
    nvec,params,type,rp,dstate,iter,elast_initial,itchat,chatter,bonddbond1,dtol,ndm)
%
% 02/06/2016
% Tim Truster
%
% Debonding interface model.
%
% type = CZM model ID
% params = [Gc sigmax dc beta ...] 
%
% elast_initial = 1 for the first iteration of a new load step, for example
% determined by transient == 0 and initia == 1; studies have shown that
% initializing the stiffness matrix with the damage tangent yields stable
% results compared to using the elastic/unloading tangent.
% dstate is a flag for what portion of the CZM model is active, e.g. damage
% or adhesion.

Ptype = 1; %TSL type
Utype = 1; %Unloading type
normsig = nvec'*tvtr;
dn = jumpu'*nvec;
tn = normsig+rp*dn;

if tn >= 0 % tension

    tvec = tvtr + rp*jumpu;

else % compression

    tvec = tvtr + rp*jumpu - tn*nvec;

end

normtvec = sqrt(tvec'*tvec);
ts = sqrt(normtvec^2 - tn^2);
tss = tvtr + rp*jumpu - tn*nvec;
if abs(ts) < 1e-12
    svec = zeros(ndm,1);
else
svec = tss/ts;
end

if kappa == 0 % bonded

 if iter > itchat && chatter == 1 % force out of chattering
    
  if tn >= 0 % tension
    
    % check yield criteria
    phi_star = phi_star_func(tn,type,params);
    
    if ts <= 1.5*phi_star % stay bonded
        
        zeta = zeros(ndm,1);
        dzdt = zeros(ndm,ndm);
        kappa = 0;
        dstate = 1;
        
    else % start debonding
        
        eta  = tvec/normtvec;
        nn = eta'*nvec;
        nt = sqrt(eta'*eta - nn^2);
        phi = phi_func(nn,nt,type,params);
        rho = rho_func(phi,0,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm);
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        
    end
    
  else % compression
      
    phi_star = phi_star_func(-1,type,params);
    
    if ts <= 1.5*phi_star % stay bonded
        
        zeta = zeros(ndm,1);
        dzdt = zeros(ndm,ndm);
        kappa = 0;
        dstate = 1;
        
    else % start debonding
        
        eta  = tvec/normtvec;
        nn = 0;
        nt = sqrt(eta'*eta - nn^2);
        phi = phi_func(nn,nt,type,params);
        rho = rho_func(phi,0,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm)-nvec*nvec';
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        
    end
      
  end
  
 else % before chattering
    
  if tn >= 0 % tension
    
    % check yield criteria
    phi_star = phi_star_func(tn,type,params);
    
    if ts <= phi_star % stay bonded
        
        zeta = zeros(ndm,1);
        dzdt = zeros(ndm,ndm);
        kappa = 0;
        bonddbond2 = 0;
        dstate = 1;
        
    else % start debonding
        
        eta  = tvec/normtvec;
        nn = eta'*nvec;
        nt = sqrt(eta'*eta - nn^2);
        phi = phi_func(nn,nt,type,params);
        rho = rho_func(phi,0,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm);
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        bonddbond2 = -1;
        
    end
    
  else % compression
      
    phi_star = phi_star_func(-1,type,params);
    
    if ts <= phi_star % stay bonded
        
        zeta = zeros(ndm,1);
        dzdt = zeros(ndm,ndm);
        kappa = 0;
        dstate = 1;
        bonddbond2 = 0;
        
    else % start debonding
        
        eta  = tvec/normtvec;
        nn = 0;
        nt = sqrt(eta'*eta - nn^2);
        phi = phi_func(nn,nt,type,params);
        rho = rho_func(phi,0,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm)-nvec*nvec';
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        bonddbond2 = -1;
        
    end
      
  end
     
 end
    
else % damaged before
    
  if tn >= 0 % tension
    
    eta  = tvec/normtvec;
    nn = eta'*nvec;
    nt = sqrt(eta'*eta - nn^2);
    phi = phi_func(nn,nt,type,params);
    dPsi = Psi_dam(kappa,Ptype,params);
    
    if (normtvec/phi - rp*kappa/phi^2 <= dPsi && ~(dstate == 4 && elast_initial)) || ...
       (normtvec/phi - rp*kappa/phi^2 <= dPsi - dtol && (dstate == 4 && elast_initial))
       % make it harder to unload during initialization step
        
        %unloading
        rho = rhoU_func(phi,kappa,normtvec,rp,Ptype,dPsi);
        zeta = rho*eta;
        Imat = eye(ndm);
        dzdt = dzdtU_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Utype,type,dPsi,params);
        dstate = 3;
        bonddbond2 = -1;
        
    else
        
        %damage
        rho = rho_func(phi,kappa,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm);
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        bonddbond2 = -1;
        
    end
    
  else % compressive
    
    eta  = tvec/normtvec;
    nn = 0;
    nt = sqrt(eta'*eta - nn^2);
    phi = phi_func(nn,nt,type,params);
    dPsi = Psi_dam(kappa,Ptype,params);
    
    if (normtvec/phi - rp*kappa/phi^2 <= dPsi && ~(dstate == 4 && elast_initial)) || ...
       (normtvec/phi - rp*kappa/phi^2 <= dPsi - dtol && (dstate == 4 && elast_initial))
       % make it harder to unload during initialization step
        
        %unloading
        rho = rhoU_func(phi,kappa,normtvec,rp,Ptype,dPsi);
        zeta = rho*eta;
        Imat = eye(ndm)-nvec*nvec';
        dzdt = dzdtU_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Utype,type,dPsi,params);
        dstate = 3;
        bonddbond2 = -1;
        
    else
        
        %damage
        rho = rho_func(phi,kappa,normtvec,rp,Ptype,params);
        zeta = rho*eta;
        kappa = rho*phi;
        Imat = eye(ndm)-nvec*nvec';
        dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params);
        dstate = 4;
        bonddbond2 = -1;
        
    end
      
  end
    
end

tvec(1) = tvec(1) - rp*zeta(1);
tvec(2) = tvec(2) - rp*zeta(2);
if (ndm == 3) % then
tvec(3) = tvec(3) - rp*zeta(3);
end

%% Determine flags for debonding
if kappa > 0 && tn >= 0
    if kappa > params(1)
        dbond = 2;
    else
        dbond = 1;
    end
elseif kappa > 0 && tn < 0
    if kappa > params(1)
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

end

function phi = phi_func(alpha,beta,type,params)
    % function for shape of yield surface
    if type == 1 % radial
        sc = params(2);
        phi = sc*sqrt(alpha^2 + beta^2);
    elseif type == 2 % ellipse
        sc = params(2);
        tc_by_sc = params(4);
        phi = sc*sqrt(alpha^2 + tc_by_sc^2*beta^2);
    end
end

function phi_star = phi_star_func(Sigma,type,params)
    % function for shape of yield surface
    if type == 1 % radial
        sc = params(2);
        if Sigma <= 0
            phi_star = sc;
        elseif Sigma <= sc
            phi_star = sqrt(sc^2 - Sigma^2);
        else
            phi_star = -1;
        end
    elseif type == 2 % ellipse
        sc = params(2);
        tc_by_sc = params(4);
        if Sigma <= 0
            phi_star = tc_by_sc*sc;
        elseif Sigma <= sc
            phi_star = tc_by_sc*sqrt(sc^2 - Sigma^2);
        else
            phi_star = -1;
        end
    end
end

function [dPsi,ddPsi] = Psi_dam(phi,Ptype,params)
    % function for type of damaging TSL
    if Ptype == 1
        Gc = params(1);
        if phi <= 2*Gc
            Psi = phi/2*(2-phi/(2*Gc));
            dPsi = 1 - phi/(2*Gc);
            ddPsi = -1/(2*Gc);
        else
            Psi = Gc;
            dPsi = 0;
            ddPsi = 0;
        end
    end
end

% function [dPsi,ddPsi] = Psi_ur(phi,Ptype,paramsur)
%     % function for elastic unloading portion
%     if Ptype == 1
%         kappa = paramsur(1);
%         dkappa = paramsur(2);
%         Psi = 1/(2*kappa)*dkappa*phi^2;
%         dPsi = phi/kappa*dkappa;
%         ddPsi = 1/kappa*dkappa;
%     end
% end

function rho = rho_func(phi,kappa,normtvec,rp,Ptype,params)
    % solve for radius of damage surface
    % rp*rho + phi^2*Psi_dam(rho,type,params) - normtvec*phi = 0
    if Ptype == 1
        Gc = params(1);
        rho = (phi^2 - normtvec*phi)/(phi^2/(2*Gc)-rp);
        rho = rho/phi;
        kappa2 = rho*phi;
        if kappa2 > 2*Gc || kappa > 2*Gc % fracture
            rho = normtvec/rp;
        end
    end
end

function dzdt = dzdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Ptype,type,params)
    % tangent matrix
    Gc = params(1);
    if kappa >= 2*Gc
        dzdt = 1/rp*Imat;
    else
        deta_dt = 1/normtvec*(Imat - eta*eta');
        drdt = drdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,nvec,svec,tn,ts,Ptype,type,params);
        dzdt = eta*drdt' + rho*deta_dt;
    end
end

function drdt = drdt_func(phi,kappa,normtvec,rp,rho,eta,ndm,nvec,svec,tn,ts,Ptype,type,params)
    % tangent matrix
    Gc = params(1);
    if Ptype == 1
        if type == 1
            dphidt = zeros(ndm,1);
        elseif type == 2
            dpda = params(2)/sqrt(tn^2+params(4)^2*ts^2)*abs(tn);
            dpdb = params(4)*params(2)/sqrt(tn^2+params(4)^2*ts^2)*abs(ts);
            phi2 = phi_func(tn,ts,type,params);
            dphidt = 1/normtvec*(dpda*nvec+dpdb*svec) - 1/normtvec^2*phi2*eta;
%             dphidt = zeros(ndm,1);
        end
        drdt = 1/(rp-phi^2/(2*Gc))*(eta-dphidt) - rho/(rp-phi^2/(2*Gc))*phi/Gc*dphidt;
    end
end

function rho = rhoU_func(phi,kappa,normtvec,rp,Utype,dPsi)
    % solve for radius of unloading surface
    if Utype == 1
        rho = normtvec/(rp+phi^2/kappa*dPsi);
    end
end

function dzdt = dzdtU_func(phi,kappa,normtvec,rp,rho,eta,ndm,Imat,nvec,svec,tn,ts,Utype,type,dPsi,params)
    % tangent matrix
    deta_dt = 1/normtvec*(Imat - eta*eta');
    drdt = drdtU_func(phi,kappa,normtvec,rp,rho,eta,ndm,nvec,svec,tn,ts,Utype,type,dPsi,params);
    dzdt = eta*drdt' + rho*deta_dt;
end

function drdt = drdtU_func(phi,kappa,normtvec,rp,rho,eta,ndm,nvec,svec,tn,ts,Utype,type,dPsi,params)
    % tangent matrix
    if Utype == 1
        if type == 1
            dphidt = zeros(ndm,1);
        elseif type == 2
            dpda = params(2)/sqrt(tn^2+params(4)^2*ts^2)*abs(tn);
            dpdb = params(4)*params(2)/sqrt(tn^2+params(4)^2*ts^2)*abs(ts);
            phi2 = phi_func(tn,ts,type,params);
            dphidt = 1/normtvec*(dpda*nvec+dpdb*svec) - 1/normtvec^2*phi2*eta;
        end
        drdt = 1/(rp+phi^2/kappa*dPsi)*eta - normtvec/(rp+phi^2/kappa*dPsi)^2*(2*phi/kappa*dPsi)*dphidt;
    end
end
