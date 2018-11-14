function [Nmat,DerN] = shf1d_Elina_M(r,ndm,ndf,p,lint)
% function [Nmat,Bmat,DetJ,c] = shf1d_Elina_M(lint,xl,lint,ndm,ndf,p,nen)
%Elina Geut 11/3/2018
%   Function is used for assembling shape functions and it's derivatives in
%   terms of glabal coordinated by multiplying 

% lint - number of nodes on current element to be considered
% ndm - number of dimensions
% lint - number of integration points used
% r - parametric variable tat gets assigned based on lint
% w - weight for numericalintegration, a;so gets assigned based on lint
% p - polynomial order
% Nmat - array of shape functions based on p

if p == 0
    for i = 1:lint
        Nmat(i,1) = 1;
        DerN(i,1) = 0;
    end      
elseif p == 1
     Nmat = zeros(ndm,ndf*lint);
     DerN = zeros(ndm,ndf*lint);
    for i=1:lint
        Nmat(i,1) = (1-r(i))/2;
        Nmat(i,2) = (1+r(i))/2;
        DerN(i,1) = -1/2;
        DerN(i,2) =1/2;
    end
elseif p == 2
     Nmat = zeros(ndm,ndf*lint);
     DerN = zeros(ndm,ndf*lint);
     
    for i=1:lint
        Nmat(i,1) = (1-r(i))/2*(-r(i));
        Nmat(i,2) =  (1+r(i))*(1-r(i));
        Nmat(i,3) = (1+r(i))/2*(r(i));
        
        DerN(i,1) = -1/2+r(i);
        DerN(i,2) = -2*r(i);
        DerN(i,3) = 1/2+r(i);
        
    end
    
elseif p == 3
        Nmat = zeros(ndm,ndf*lint);
        DerN = zeros(ndm,ndf*lint);
    for i=1:lint
        Nmat(i,1) =  1/16*(9*r(i)^2 - 1)*(1 - r(i));
        Nmat(i,2) =  9/16*(r(i)^2 - 1)*(3*r(i) - 1);
        Nmat(i,3) =  9/16*(1 - r(i)^2)*(3*r(i) + 1);
        Nmat(i,4) =  1/16*(9*r(i)^2 - 1)*(r(i) + 1);
            
        DerN(i,1) =  1/16*(1 + r(i)*(18 - 27*r(i)));
        DerN(i,2) =  9/16*(-3 - r(i)*(2 - 9*r(i)));
        DerN(i,3) =  9/16*(3 - r(i)*(2 + 9*r(i)));
        DerN(i,4) =  1/16*(-1 + r(i)*(18 + 27*r(i)));

    end
         
end