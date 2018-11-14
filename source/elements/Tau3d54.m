function [tau, hr1, intb, Cmat] = Tau3d54(xlint,xl,ul,ul_n,hr,mateprop,plasmodel,plasversion,nel,nen,ndf,forceelast)
% 05/10/2015
% computes tau for edge bubble for von Mises plasticity, NL_Elem53_3d
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
One = [1; 1; 1; 0; 0; 0];

% forceelast = 1; % set true to force elastic moduli for tau

tau = 0;
ib = 0;
der = 0;
bf = 1;
intb = 0;
    if nel == 4
        lint = 4;11;36;
    elseif nel == 8
        lint = 8;1000;64;27;
    elseif nel == 10
        lint = 24;36;
    else
        lint = 27;1000;
    end
Bmat = zeros(7,nel*ndf);
hr1 = zeros(13,1);

    
% Evaluate plastic variables at center of ELEMENT
if nel == 4 || nel == 10
[w,ss] =  int3d_t(1,1,ib);
[shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
QXY = shgtt(xl,nel,shld,shls,nen,bf,der,be);
elseif nel == 8 || nel == 27
[w,ss] = intpntb(1,1,ib);
[shl,shld,shls,be] = shlq(ss,nel,nel,der,bf);
QXY = shgq(xl,nel,shld,shls,nen,bf,der,be); 
end

% Form B matrix for coarse-scale strains
for ie = 1:nel
Bmat(:,(ie-1)*4+1:4*ie) = [QXY(ie,1) 0        0 0
                           0        QXY(ie,2) 0 0
                           0        0 QXY(ie,3) 0
                           QXY(ie,2) QXY(ie,1) 0 0
                           0 QXY(ie,3) QXY(ie,2) 0
                           QXY(ie,3) 0 QXY(ie,1) 0
                           0 0 0 0];
end

if forceelast
    
    % Compute input for Radial Return from coarse-scale fields
    du = 0*Bmat*reshape(ul,nel*ndf,1);
    eps3d = du(1:6); %2D enhanced strain
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 6;
    ahr = 13;
    ep_n = 0*[hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
    beta_n = 0*[hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
    a_n = 0*hr(ahr);

    [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

    % Convert output from 3D to 2D
    Cmat = Cdev_n1;
    
else
    
    % Compute input for Radial Return from coarse-scale fields
    du = 0*Bmat*reshape(ul,nel*ndf,1);
    eps3d = du(1:6); %2D enhanced strain
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 6;
    ahr = 13;
    ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
    beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
    a_n = hr(ahr);

    [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0Tau(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

    % Convert output from 3D to 2D
    Cmat = Cdev_n1;

    % Store history variables
    ephr = 0; %pointer for plastic strain at pt l
    betahr = 6;
    ahr = 13;
    hr1(ephr+1:ephr+6) = ep_n1(1:6);
    hr1(betahr+1:betahr+6) = beta_n1(1:6);
    hr1(ahr) = a_n1; 

end % forceelast


% Loop to compute bubble function energy integral
for ll = 1:lint
    
    % Evaluate integ weight, determinant, and Jacobian over truncated
    % sector
    if nel == 4 || nel == 10
      [w,ss] =  int3d_t(ll,lint,ib);
      [~,shld,shls,be] = shltt(ss,nel,nel,der,bf);
      [~,~,Jdet,~,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);

    else
      [w,ss] =  intpntb(ll,lint,ib);
      [~,shld,shls,be] = shlb(ss,nel,nel,der,bf);
      [~,~,Jdet,~,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);

    end
    
    % Evaluate bubble function over master sector
    if nel == 4 || nel == 10
    [b,db] = facebubbleT(ss,nel); 
    else
    [b,db] = facebubbleQ(ss,nel);   
    end
    %Map derivatives forward to truncated sector reference configuration
    bee = (db'/xs);
         
	xsj = Jdet;   
        
    Bmatb = [bee(1)  0       0
            0       bee(2)  0
            0       0       bee(3)
            bee(2)  bee(1)  0
            0       bee(3)  bee(2)
            bee(3)  0       bee(1)];
            
    tau = tau + w*xsj*(Bmatb'*Cmat*Bmatb);

    intb = intb + xsj*w*b;


end

tau = inv(tau);
