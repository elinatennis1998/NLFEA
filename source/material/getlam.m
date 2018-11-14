function lam = getlam(mateprop)
%
% Compute lamda for volumetric component of material model
% 10/16/2011

mati = mateprop(1);

if mati <= 4
    if mati <= 2
        if mati == 1 % Standard Neo-Hookean Material

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));%4e5;
            
        else % mati == 2 % Isochoric Neo-Hookean Material

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = PatchE/(3*(1-2*Patchv));%2.8333e3;%501;%40.0942e4;%

        end
    elseif mati == 3 % Mooney-Rivlin Material

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        
    else % mati == 4 % Ogden Material, Reese IJNME38
        
            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
            
    end
elseif mati > 5 && mati <= 8
    if mati <= 6
        if mati == 5 % Yeoh Material, Stein CMAME190
        
            PatchE = mateprop(3);
            Patchv = mateprop(4);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
            
        else % mati == 6
        
            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
            
        end
    elseif mati == 7 % Ogden Material, Simo CMAME85
        
        PatchE = mateprop(4);
        Patchv = mateprop(5);
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        
    else % mati == 8 % Hencky model

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = PatchE/(3*(1-2*Patchv));
        
    end
    
else
        if mati == 9 % Compressible anisotropic

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));%4e5;
            
        elseif mati == 10 % Compressible double anisotropic

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));%4e5;
            
        end 
end