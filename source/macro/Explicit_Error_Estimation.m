% 06/15/2013
% Tim Truster
% Explicit Error Calculation, segregated from NL_FEA_Program
            
if numEn == 16
    
    bubblevals = zeros(numel,2);
    Ieffvals = zeros(numel,5);
    
    darcystokeserror = 0;
    
    isw = 11;
    FormFE
    
    Err = Energy(1:9);
    Fine = Energy(10:15);
%     Fine = Energy(10:16);
    L2uerr = sqrt(Err(1) + Err(2));
    L2uerr = log(L2uerr)/log(10);
    L2perr = sqrt(Err(3));
    L2perr = log(L2perr)/log(10);
    H1uerr = sqrt(Err(4) + Err(5) + Err(7) + Err(8));
    H1uerr = log(H1uerr)/log(10);
    H1perr = sqrt(Err(6) + Err(9));
    H1perr = log(H1perr)/log(10);
    L2Fine = log(sqrt(Fine(1)+Fine(2)))/log(10);
%     L2pFine = log(sqrt(Fine(7)))/log(10);
    H1Fine = log(sqrt(Fine(3)+Fine(4)+Fine(5)+Fine(6)))/log(10);

%     h = log(1)/log(10);
    fprintf('Standard %i  %1.7e  %1.7e  %1.7e  %1.7e  %i\n',numel,L2uerr,L2perr,H1uerr,H1perr,iprob)
    fprintf('Explicit %i  %1.7e  %1.7e\n',numel,L2Fine,H1Fine)
    
elseif numEn == 6
    
    bubblevals = zeros(numel,2);
    Ieffvals = zeros(numel,5);
    
    darcystokeserror = 0;
    
    isw = 11;
    FormFE
    
    Err = Energy(1:6);
    L2uerr = sqrt(Err(1) + Err(2));
    L2uerr = log(L2uerr)/log(10);
    H1uerr = sqrt(Err(3) + Err(4) + Err(5) + Err(6));
    H1uerr = log(H1uerr)/log(10);

%     h = log(1)/log(10);
    fprintf('Standard %i  %1.7e  %1.7e  %i\n',numel,L2uerr,H1uerr,iprob)
    
elseif numEn == 30
    
    bubblevals = zeros(numel,3);
    Ieffvals = zeros(numel,5);
    
    darcystokeserror = 0;
    
    isw = 11;
    FormFE
    
    Err = Energy(1:16);
    L2uerr = sqrt(Err(1) + Err(2) + Err(3));
    L2uerr = log(L2uerr)/log(10);
    L2perr = sqrt(Err(4));
    L2perr = log(L2perr)/log(10);
    H1uerr = sqrt(Err(5) + Err(6) + Err(7) + Err(9) + Err(10) + Err(11) + Err(13) + Err(14) + Err(15));
    H1uerr = log(H1uerr)/log(10);
    H1perr = sqrt(Err(8) + Err(12) + Err(16));
    H1perr = log(H1perr)/log(10);

%     h = log(1)/log(10);
    fprintf('Standard %i  %1.7e  %1.7e  %1.7e  %1.7e  %i\n',numel,L2uerr,L2perr,H1uerr,H1perr,iprob)
    
elseif numEn == 12

    Ieffvals = zeros(numel,5);
    
    isw = 11;
    FormFE
    
    Err = Energy(1:12);
    L2uerr = sqrt(Err(1) + Err(2) + Err(3));
    L2uerr = log(L2uerr)/log(10);
    H1uerr = sqrt(Err(4) + Err(5) + Err(6) + Err(7) + Err(8) + Err(9) + Err(10) + Err(11) + Err(12));
    H1uerr = log(H1uerr)/log(10);

%     h = log(1)/log(10);
    fprintf('Standard %i  %1.7e  %1.7e  %i\n',numel,L2uerr,H1uerr,iprob)
    
elseif numEn == 2
    
    Ieffvals = zeros(numel,5);
    
    isw = 11;
    FormFE
    
    L2uerr = sqrt(Energy(1));
    L2uerr = log(L2uerr)/log(10);
    H1uerr = sqrt(Energy(2));
    H1uerr = log(H1uerr)/log(10);
    fprintf('Standard %i  %1.7e  %1.7e  %i\n',numel,L2uerr,H1uerr,iprob)
         
else

    Ieffvals = zeros(numel,5);
    
    isw = 11;
    FormFE

    Err = Energy
        
end