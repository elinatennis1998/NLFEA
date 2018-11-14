%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke,Re,Ie]=brick_3d(e_mat,e_spa,i_var,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stiffness matrix Ke and righthand side Re for quads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nmat,ngrw] = size(mat);
 nod  = 8;
 ndim = 3;
 nip  = 2;

Ie = i_var;
Re = zeros(nod*ndim,1); 
Ke = zeros(nod*ndim,nod*ndim); 

indx=[1;4;7;10;13;16;19;22]; ex_mat=e_mat(indx); 
indy=[2;5;8;11;14;17;20;23]; ey_mat=e_mat(indy); 
indz=[3;6;9;12;15;18;21;24]; ez_mat=e_mat(indz); 

JT =zeros(ndim*nip*nip*nip,ndim);
N  =zeros(nip*nip*nip,nod);
dNr=zeros(ndim*nip*nip*nip,nod);
wp =zeros(nip*nip*nip,1);

% integration points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1=0.577350269189626; w1=1;
gp(:,1)=[-g1; g1;-g1; g1;-g1; g1;-g1; g1];  
gp(:,2)=[-g1;-g1; g1; g1;-g1;-g1; g1; g1];
gp(:,3)=[-g1;-g1;-g1;-g1; g1; g1; g1; g1];
 w(:,1)=[ w1; w1; w1; w1; w1; w1; w1; w1];   
 w(:,2)=[ w1; w1; w1; w1; w1; w1; w1; w1];
 w(:,3)=[ w1; w1; w1; w1; w1; w1; w1; w1];

wp   = w(:,1).*w(:,2).*w(:,3); 
xsi = gp(:,1);  
eta = gp(:,2);  
rho = gp(:,3);
r3   = ndim*nip*nip*nip;

% shape functions @ 8 integration points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
N(:,1)=(1-xsi).*(1-eta).*(1-rho)/8;  N(:,2)=(1+xsi).*(1-eta).*(1-rho)/8;
N(:,3)=(1+xsi).*(1+eta).*(1-rho)/8;  N(:,4)=(1-xsi).*(1+eta).*(1-rho)/8;
N(:,5)=(1-xsi).*(1-eta).*(1+rho)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+rho)/8;
N(:,7)=(1+xsi).*(1+eta).*(1+rho)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+rho)/8;

% partial derivatives of shape functions wrt xsi & eta & rho %%%%%%%%%%%
dNr(1:3:r3,1)=-(1-eta).*(1-rho)/8; dNr(1:3:r3,2)=+(1-eta).*(1-rho)/8;
dNr(1:3:r3,3)=+(1+eta).*(1-rho)/8; dNr(1:3:r3,4)=-(1+eta).*(1-rho)/8;
dNr(1:3:r3,5)=-(1-eta).*(1+rho)/8; dNr(1:3:r3,6)=+(1-eta).*(1+rho)/8;
dNr(1:3:r3,7)=+(1+eta).*(1+rho)/8; dNr(1:3:r3,8)=-(1+eta).*(1+rho)/8;

dNr(2:3:r3,1)=-(1-xsi).*(1-rho)/8; dNr(2:3:r3,2)=-(1+xsi).*(1-rho)/8;
dNr(2:3:r3,3)=+(1+xsi).*(1-rho)/8; dNr(2:3:r3,4)=+(1-xsi).*(1-rho)/8;
dNr(2:3:r3,5)=-(1-xsi).*(1+rho)/8; dNr(2:3:r3,6)=-(1+xsi).*(1+rho)/8;
dNr(2:3:r3,7)=+(1+xsi).*(1+rho)/8; dNr(2:3:r3,8)=+(1-xsi).*(1+rho)/8;

dNr(3:3:r3,1)=-(1-xsi).*(1-eta)/8; dNr(3:3:r3,2)=-(1+xsi).*(1-eta)/8;
dNr(3:3:r3,3)=-(1+xsi).*(1+eta)/8; dNr(3:3:r3,4)=-(1-xsi).*(1+eta)/8;
dNr(3:3:r3,5)=+(1-xsi).*(1-eta)/8; dNr(3:3:r3,6)=+(1+xsi).*(1-eta)/8;
dNr(3:3:r3,7)=+(1+xsi).*(1+eta)/8; dNr(3:3:r3,8)=+(1-xsi).*(1+eta)/8;

JT=dNr*[ex_mat;ey_mat;ez_mat]';

for ip=1:(nip*nip*nip)
  indx=[ndim*ip-2; ndim*ip-1; ndim*ip];
  if (det(JT(indx,:))<10*eps); disp('det(J) equal or less than zero!'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop over all integration points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ip=1:(nip*nip*nip)
  indx=[3*ip-2; 3*ip-1; 3*ip];
  detJ  = det(JT(indx,:));
  JTinv = inv(JT(indx,:));
  dNx   = JTinv*dNr(indx,:);

  F=zeros(ndim,ndim);
  for j=1:nod
    jndx=[3*j-2; 3*j-1; 3*j];
    F=F+e_spa(jndx)'*dNx(:,j)';
  end
  
  var = i_var(ip);  
  if (ngrw==7); [A,P,var]=cnst_den3(F,var,mat,ndim); end;  % density
  if (ngrw==9); [A,P,var]=cnst_vol(F,var,mat,ndim); end;  % volume
  Ie(ip)=var;
    
%%% righthand sides -- loop over all nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nod
  eni=(i-1)*ndim;
  Re(eni+1)=Re(eni+1)+(P(1,1)*dNx(1,i)'+P(1,2)*dNx(2,i)'+P(1,3)*dNx(3,i)')*detJ*wp(ip);
  Re(eni+2)=Re(eni+2)+(P(2,1)*dNx(1,i)'+P(2,2)*dNx(2,i)'+P(2,3)*dNx(3,i)')*detJ*wp(ip);   
  Re(eni+3)=Re(eni+3)+(P(3,1)*dNx(1,i)'+P(3,2)*dNx(2,i)'+P(3,3)*dNx(3,i)')*detJ*wp(ip); 
end  

%%% stiffness matrices -- loop over all nodes %%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nod
for j=1:nod
  eni=(i-1)*ndim;
  enj=(j-1)*ndim;

  Ke(eni+1,enj+1)=Ke(eni+1,enj+1) ...
                 +(dNx(1,i)*(A(1,1,1,1)*dNx(1,j)+A(1,1,1,2)*dNx(2,j)+A(1,1,1,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(1,2,1,1)*dNx(1,j)+A(1,2,1,2)*dNx(2,j)+A(1,2,1,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(1,3,1,1)*dNx(1,j)+A(1,3,1,2)*dNx(2,j)+A(1,3,1,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+1,enj+2)=Ke(eni+1,enj+2) ...
                 +(dNx(1,i)*(A(1,1,2,1)*dNx(1,j)+A(1,1,2,2)*dNx(2,j)+A(1,1,2,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(1,2,2,1)*dNx(1,j)+A(1,2,2,2)*dNx(2,j)+A(1,2,2,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(1,3,2,1)*dNx(1,j)+A(1,3,2,2)*dNx(2,j)+A(1,3,2,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+1,enj+3)=Ke(eni+1,enj+3) ...
                 +(dNx(1,i)*(A(1,1,3,1)*dNx(1,j)+A(1,1,3,2)*dNx(2,j)+A(1,1,3,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(1,2,3,1)*dNx(1,j)+A(1,2,3,2)*dNx(2,j)+A(1,2,3,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(1,3,3,1)*dNx(1,j)+A(1,3,3,2)*dNx(2,j)+A(1,3,3,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+2,enj+1)=Ke(eni+2,enj+1) ...
                 +(dNx(1,i)*(A(2,1,1,1)*dNx(1,j)+A(2,1,1,2)*dNx(2,j)+A(2,1,1,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(2,2,1,1)*dNx(1,j)+A(2,2,1,2)*dNx(2,j)+A(2,2,1,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(2,3,1,1)*dNx(1,j)+A(2,3,1,2)*dNx(2,j)+A(2,3,1,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+2,enj+2)=Ke(eni+2,enj+2) ...
                 +(dNx(1,i)*(A(2,1,2,1)*dNx(1,j)+A(2,1,2,2)*dNx(2,j)+A(2,1,2,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(2,2,2,1)*dNx(1,j)+A(2,2,2,2)*dNx(2,j)+A(2,2,2,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(2,3,2,1)*dNx(1,j)+A(2,3,2,2)*dNx(2,j)+A(2,3,2,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+2,enj+3)=Ke(eni+2,enj+3) ...
                 +(dNx(1,i)*(A(2,1,3,1)*dNx(1,j)+A(2,1,3,2)*dNx(2,j)+A(2,1,3,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(2,2,3,1)*dNx(1,j)+A(2,2,3,2)*dNx(2,j)+A(2,2,3,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(2,3,3,1)*dNx(1,j)+A(2,3,3,2)*dNx(2,j)+A(2,3,3,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+3,enj+1)=Ke(eni+3,enj+1) ...
                 +(dNx(1,i)*(A(3,1,1,1)*dNx(1,j)+A(3,1,1,2)*dNx(2,j)+A(3,1,1,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(3,2,1,1)*dNx(1,j)+A(3,2,1,2)*dNx(2,j)+A(3,2,1,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(3,3,1,1)*dNx(1,j)+A(3,3,1,2)*dNx(2,j)+A(3,3,1,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+3,enj+2)=Ke(eni+3,enj+2) ...
                 +(dNx(1,i)*(A(3,1,2,1)*dNx(1,j)+A(3,1,2,2)*dNx(2,j)+A(3,1,2,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(3,2,2,1)*dNx(1,j)+A(3,2,2,2)*dNx(2,j)+A(3,2,2,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(3,3,2,1)*dNx(1,j)+A(3,3,2,2)*dNx(2,j)+A(3,3,2,3)*dNx(3,j)))*detJ*wp(ip);
  Ke(eni+3,enj+3)=Ke(eni+3,enj+3) ...
                 +(dNx(1,i)*(A(3,1,3,1)*dNx(1,j)+A(3,1,3,2)*dNx(2,j)+A(3,1,3,3)*dNx(3,j)) ...
                 + dNx(2,i)*(A(3,2,3,1)*dNx(1,j)+A(3,2,3,2)*dNx(2,j)+A(3,2,3,3)*dNx(3,j)) ...
                 + dNx(3,i)*(A(3,3,3,1)*dNx(1,j)+A(3,3,3,2)*dNx(2,j)+A(3,3,3,3)*dNx(3,j)))*detJ*wp(ip);     
                        
end, end

end
% loop over all integration points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
