%********************************************************************************************* 
function [F,dF]=dexpmt(M,dM) % by Lubomír Bran?ík, 2008 
% Scale M by power of 2 so that its norm is < 1/2 
[f,e]=log2(norm(M,'inf')); 
r=max(0,e+1); 
M=M/2^r; dM=dM/2^r; 
% Taylor series expansion of exp(M) and diff[exp(M)] 
FM=eye(size(M)); dFM=dM; 
F=FM+M; dF=dM; 
l=1; n=1; 
incF=1; 
while incF>1e-16 
 n=n+1; l=l*n; 
 FM=FM*M; 
 dFM=dM*FM+M*dFM; 
 dF=dF+dFM/l; 
 F=F+FM*M/l; 
 incF=sum(sum(abs(FM)+abs(dFM))); 
end 
% Undo scaling by repeated squaring 
for k=1:r 
 dF=dF*F+F*dF; 
 F=F*F; 
end 
%*********************************************************************************************