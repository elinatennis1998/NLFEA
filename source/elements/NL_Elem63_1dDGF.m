% 11/4/2016


%% Task Switch - Perform various element functions
switch isw 
    case 1
nh1= 2*1; % number of variables=2, number of integration points=1
%%
case 3 % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = NL_Elem63_1dDGF03(mateprop,ul,xl,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
end