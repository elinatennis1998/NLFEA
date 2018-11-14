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
ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c1*(Nmat1'*IntFor);
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
ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c2*(Nmat2'*IntFor);
ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) = ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) ...
    + c2*(JxXif*Nmat2'*kIF/tstep*Nmat2);
