function [Bmat,DetJ] = numint_Elina(xl,DerN,nel)
%Elia Geut. Last modified 11/7/2018
%   This function performs numrical integration 
%Change the shape functions
% J - Jacobain of the element 
% Bmat - derivatives of Nmat
% xl - biunit coordinates of a current element 
% nel - number of nodes on element
 J = DerN*xl(:,1:nel)';
 Bmat = (1/J)*DerN;
 DetJ = abs(det(J));
 
end



