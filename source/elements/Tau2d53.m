function [tau, hr1, intb] = Tau2d53(xlint,xl,ul,ul_n,hr,mateprop,plasmodel,plasversion,nel,nen,ndf,forceelast)
% 05/12/2014
% computes tau for edge bubble for von Mises plasticity, NL_Elem53_2d
%
% uses proper treatment of truncated sectors
%
% Input is history and displacement from coarse scales
% Output is stability tensor and updated history for coarse scales

% Set Material Properties

ElemYM = mateprop(2);
Elemv = mateprop(3);
Khard = mateprop(4);
Hhard = mateprop(5);
sigy = mateprop(6);
bulk = ElemYM/(3*(1-2*Elemv));
mu = ElemYM/(2*(1+Elemv));
One = [1; 1; 0];
ind3to2 = [1 2 4];

% forceelast = 1;

tau = 0;
ib = 0;
der = 0;
bf = 1;
intb = 0;
if nel == 3 % Matches with Table 3 in the NLDG paper
    lint = 4;
elseif nel == 4
    lint = 9;
elseif nel == 6
    lint = 13;
else
    lint = 16;
end
Bmat = zeros(3,nel*ndf);
hr1 = zeros(7,1);

    
% Evaluate plastic variables at center of ELEMENT
if nel == 3 || nel == 6
[w,litr,lits] = intpntt(1,1,ib);
[shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
QXY = shgt(xl,nel,shld,shls,nen,bf,der,be);
elseif nel == 4 || nel == 9
[w,litr,lits] = intpntq(1,1,ib);
[shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
QXY = shgq(xl,nel,shld,shls,nen,bf,der,be); 
end
    
%     % Find physical location of integration point
%     xy = xl*shl;
%     POUxi = POU_Coord(xy(1),xy(2),xlint,0,nel);
%     % Evaluate bubble function at physical location in interface SECTOR
%     if nel == 3 || nel == 6
%     [shl,shld,shls,be] = shlt(POUxi(1),POUxi(2),3,3,der,bf);
%     [~, ~, ~, ~, xs] = shgt(xlint,3,shld,shls,3,bf,der,be);
%     elseif nel == 4 || nel == 9
%     [shl,shld,shls,be] = shlq(POUxi(1),POUxi(2),4,4,der,bf);
%     [~, ~, ~, ~, xs] = shgq(xl,4,shld,shls,4,bf,der,be); 
%     end
%     
%     % Evaluate bubble function over master sector
%     if nel == 3 || nel == 6
%     [~,db] = edgebubble(litr,lits,nel);
%     else
%     [~,db] = edgebubbleQ(litr,lits,nel);
%     end
%     %Map derivatives forward to truncated sector reference configuration
%     db = (db'/xs);

% Form B matrix for coarse-scale strains
for ie = 1:nel
Bmat(:,(ie-1)*2+1:2*ie) = [QXY(ie,1) 0       
                           0        QXY(ie,2)
                           QXY(ie,2) QXY(ie,1)];
end

if forceelast
    
    % Compute input for Radial Return from coarse-scale fields
    du = 0*Bmat*reshape(ul,nel*ndf,1);
    eps2d = du(1:3); %2D enhanced strain
    eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 3;
    ahr = 7;
    ep_n = 0*[hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
    ep_n(3) = -(ep_n(1) + ep_n(2));
    beta_n = 0*[hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
    beta_n(3) = -(beta_n(1) + beta_n(2));
    a_n = 0*hr(ahr);

    [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

    % Convert output from 3D to 2D
    Cmat = C_n1(ind3to2,ind3to2);
    
else

if plasversion == 1
    
    % Compute input for Radial Return from coarse-scale fields
    eps2d = Bmat*reshape(ul,nel*ndf,1); %2D strain
    eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 3;
    ahr = 7;
    ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
    ep_n(3) = -(ep_n(1) + ep_n(2));
    beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
    beta_n(3) = -(beta_n(1) + beta_n(2));
    a_n = hr(ahr);

    [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

    % Convert output from 3D to 2D
    Cmat = C_n1(ind3to2,ind3to2);

    % Store history variables
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 3;
    ahr = 7;
    hr1(ephr+1) = ep_n1(1);
    hr1(ephr+2) = ep_n1(2);
    hr1(ephr+3) = ep_n1(4);
    hr1(betahr+1) = beta_n1(1);
    hr1(betahr+2) = beta_n1(2);
    hr1(betahr+3) = beta_n1(4);
    hr1(ahr) = a_n1; 
                
else
    
    % Compute input for Radial Return
    deps2d = Bmat*reshape(ul-ul_n,nel*ndf,1); %3D enhanced strain
    deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
    ephr = 0; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
    ahr = 7;
    ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
    a_n = hr(ahr);
    eps_n1 = ee_n + deps_n1;

    if plasmodel == 1
    [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
    C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
    elseif plasmodel == 2
    [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
    C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
    end

    % Convert output from 3D to 2D
    Cmat = C_n1(1:3,1:3);

    % Store history variables
    ephr = 0; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
    ahr = 7;
    hr1(ephr+1) = RSTAVA(1);
    hr1(ephr+2) = RSTAVA(2);
    hr1(ephr+3) = RSTAVA(3);
    hr1(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
    hr1(ahr) = RSTAVA(5);
    
end

end %forceelast

% Loop to compute bubble function energy integral
for ll = 1:lint
    
    % Evaluate integ weight, determinant, and Jacobian over truncated
    % sector
    if nel == 3 || nel == 6
    [w,litr,lits] = intpntt(ll,lint,ib);
    [~,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
    [~, ~, Jdet, ~, xs] = shgt(xlint,3,shld,shls,3,bf,der,be);
    elseif nel == 4 || nel == 9
    [w,litr,lits] = intpntq(ll,lint,ib);
    [~,shld,shls,be] = shlq(litr,lits,4,4,der,bf);
    [~, ~, Jdet, ~, xs] = shgq(xlint,4,shld,shls,4,bf,der,be); 
    end
    
    % Evaluate bubble function over master sector
    if nel == 3 || nel == 6
    [b,db] = edgebubble(litr,lits,nel);
    else
    [b,db] = edgebubbleQ(litr,lits,nel);
    end
    %Map derivatives forward to truncated sector reference configuration
    db = (db'/xs);
         
	xsj = Jdet;   
     
    Bmatb = [db(1)  0       
             0      db(2)  
             db(2)  db(1)];
            
    tau = tau + w*xsj*(Bmatb'*Cmat*Bmatb);

    intb = intb + xsj*w*b;


end

tau = inv(tau);
