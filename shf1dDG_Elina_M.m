function [NmatL,DerNL,NmatR,DerNR] = shf1dDG_Elina_M(rL,rR,ndm,ndf,pL,pR,lintL,lintR)
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

if pL == 0 
    for i = 1:lintL
        NmatL(i,1) = 1;
        DerNL(i,1) = 0;
    end      
elseif pL == 1
     NmatL = zeros(ndm,ndf*lintL);
     DerNL = zeros(ndm,ndf*lintL);
    for i=1:lintL
        NmatL(i,1) = (1-rL(i))/2;
        NmatL(i,2) = (1+rL(i))/2;
        DerNL(i,1) = -1/2;
        DerNL(i,2) =1/2;
    end
elseif pL == 2
     NmatL = zeros(ndm,ndf*lintL);
     DerNL = zeros(ndm,ndf*lintL);
     
    for i=1:lintL
        NmatL(i,1) = (1-rL(i))/2*(-rL(i));
        NmatL(i,2) =  (1+rL(i))*(1-rL(i));
        NmatL(i,3) = (1+rL(i))/2*(rL(i));
        
        DerNL(i,1) = -1/2+rL(i);
        DerNL(i,2) = -2*rL(i);
        DerNL(i,3) = 1/2+rL(i);
        
    end
end
    
    
if pR == 0 
    for i = 1:lintR
        NmatR(i,1) = 1;
        DerNR(i,1) = 0;
    end      
elseif pR == 1
     NmatR = zeros(ndm,ndf*lintR);
     DerNR = zeros(ndm,ndf*lintR);
    for i=1:lintR
        NmatR(i,1) = (1-rR(i))/2;
        NmatR(i,2) = (1+rR(i))/2;
        DerNR(i,1) = -1/2;
        DerNR(i,2) =1/2;
    end
elseif pR == 2
     NmatR = zeros(ndm,ndf*lintR);
     DerNR = zeros(ndm,ndf*lintR);
     
    for i=1:lintR
        NmatR(i,1) = (1-rR(i))/2*(-rR(i));
        NmatR(i,2) =  (1+rR(i))*(1-rR(i));
        NmatR(i,3) = (1+rR(i))/2*(rR(i));
        
        DerNR(i,1) = -1/2+rR(i);
        DerNR(i,2) = -2*rR(i);
        DerNR(i,3) = 1/2+rR(i);
    end
end 
    
% elseif p == 3
%         Nmat = zeros(ndm,ndf*lint);
%         DerN = zeros(ndm,ndf*lint);
%     for i=1:lint
%         Nmat(i,1) =  1/16*(9*r(i)^2 - 1)*(1 - r(i));
%         Nmat(i,2) =  9/16*(r(i)^2 - 1)*(3*r(i) - 1);
%         Nmat(i,3) =  9/16*(1 - r(i)^2)*(3*r(i) + 1);
%         Nmat(i,4) =  1/16*(9*r(i)^2 - 1)*(r(i) + 1);
%             
%         DerN(i,1) =  1/16*(1 + r(i)*(18 - 27*r(i)));
%         DerN(i,2) =  9/16*(-3 - r(i)*(2 - 9*r(i)));
%         DerN(i,3) =  9/16*(3 - r(i)*(2 + 9*r(i)));
%         DerN(i,4) =  1/16*(-1 + r(i)*(18 + 27*r(i)));
% 
%     end
         

end

