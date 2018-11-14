% Tim Truster
% 08/03/2016
%
% Evaluate interactive force for drag type, called by NL_Elem34_3dBoth;
% written from InterForce3 and 4.

IntFor = zeros(ndf,1);

% get values of constituents from both sides
nelA = nel1;
nelB = nel2;
xl1 = xl(1:ndm,1:nel1);
ul1 = ul(1:ndf,1:nel1);
ul1_n = ul_n(1:ndf,1:nel1);
xl2 = xl(1:ndm,nel1+1:nel1+nel2);
ul2 = ul(1:ndf,nel1+1:nel1+nel2);
ul2_n = ul_n(1:ndf,nel1+1:nel1+nel2);


%% First interaction
xlA = xl1;
ulA = ul1;
ulA_n = ul1_n;
xlB = xl2;
ulB = ul2;
ulB_n = ul2_n;
JxXif = JxX1;
    
% Current physical location of int pt
xint = (xlA(1,1:nelA)+ulA(1,1:nelA))*shl;
yint = (xlA(2,1:nelA)+ulA(2,1:nelA))*shl;
zint = (xlA(3,1:nelA)+ulA(3,1:nelA))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xlB(:,1:nelB)+ulB(1:3,1:nelB),1,nelA);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nelB,nelB,0,0);
v1 = (ulA(1:3,1:nelA)-ulA_n(1:3,1:nelA))/tstep*shl;
v2 = (ulB(1:3,1:nelB)-ulB_n(1:3,1:nelB))/tstep*shl2;

IntFor = JxXif*kIF*(v1 - v2);

% Copied from the single elements NL_Elem26_3dMG and NL_Elem33_3dM
ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c1*(Nmat1'*JxXif*kIF*v1);
ElemK(1:nel1*ndf,1:nel1*ndf) = ElemK(1:nel1*ndf,1:nel1*ndf) ...
    + c1*(JxXif*Nmat1'*kIF/tstep*Nmat1);


%% Second interaction
xlA = xl2;
ulA = ul2;
ulA_n = ul2_n;
xlB = xl1;
ulB = ul1;
ulB_n = ul1_n;
JxXif = JxX2;
v2 = (ulA(1:3,1:nelA)-ulA_n(1:3,1:nelA))/tstep*shl;

ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c2*(Nmat2'*JxXif*kIF*v2);
ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) = ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) ...
    + c2*(JxXif*Nmat2'*kIF/tstep*Nmat2);


%% Projection terms
% 1-2
    
v2 = (ul2(1:3,1:nel2)-ul2_n(1:3,1:nel2))/tstep*shl;

% Current physical location of int pt
xint = (xl2(1,1:nel2)+ul2(1,1:nel2))*shl;
yint = (xl2(2,1:nel2)+ul2(2,1:nel2))*shl;
zint = (xl2(3,1:nel2)+ul2(3,1:nel2))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl1(:,1:nel1)+ul1(1:3,1:nel1),1,nel1);

% Evaluate  basis functions at integration point
shl1 = shlb(xi,nel1,nel1,0,0);
Nmat1p = zeros(3,nel*ndf/2);
for mm = 1:nel1
  Nmat1p(:,3*mm-2:3*mm) = [shl1(mm) 0        0         
                           0        shl1(mm) 0         
                           0        0        shl1(mm)];
end

ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c2*(-Nmat1p'*JxX2*kIF*v2);

ElemK(1:nel1*ndf,nel1*ndf+1:nel*ndf) = ElemK(1:nel1*ndf,nel1*ndf+1:nel*ndf) ...
    - c2*(Nmat1p'*JxX2*kIF/tstep*Nmat2);


% 2-1
    
v1 = (ul1(1:3,1:nel1)-ul1_n(1:3,1:nel1))/tstep*shl;

% Current physical location of int pt
xint = (xl1(1,1:nel1)+ul1(1,1:nel1))*shl;
yint = (xl1(2,1:nel1)+ul1(2,1:nel1))*shl;
zint = (xl1(3,1:nel1)+ul1(3,1:nel1))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl2(:,1:nel2)+ul2(1:3,1:nel2),1,nel2);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nel2,nel2,0,0);
Nmat2p = zeros(3,nel*ndf/2);
for mm = 1:nel1
  Nmat2p(:,3*mm-2:3*mm) = [shl2(mm) 0        0         
                           0        shl2(mm) 0         
                           0        0        shl2(mm)];
end

ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c1*(-Nmat2p'*JxX1*kIF*v1);

ElemK(nel1*ndf+1:nel*ndf,1:nel1*ndf) = ElemK(nel1*ndf+1:nel*ndf,1:nel1*ndf) ...
    - c1*(Nmat2p'*JxX1*kIF/tstep*Nmat1);
         