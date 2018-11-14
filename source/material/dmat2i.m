function [dmati] = dmat2i(mateprop)
%
% 6th order Shear constitutive material model evaluation for NL elasticity
% 6/17/2013
% Pinlei Chen
% Output:
%         D = F_iI*F_jJ*F_kK*F_lL*F_mM*F_nN*D_IJKLMN
% easier to go to 3D

mati = mateprop(1);

if mati <= 4
    if mati <= 2
        if mati == 1 % Standard Neo-Hookean Material
            
        dpmat3 = [8 0 0 
                  0 0 0 
                  0 0 2 
                  0 0 0 
                  0 8 0 
                  0 0 2 
                  0 0 2 
                  0 0 2 
                  2 2 0 ];

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            mu = PatchE/(2*(1+Patchv));%80.19;
            dmati=mu*dpmat3;
        end
    end
else
    if mati <= 6
        if mati == 6 % Anisotropic, Stein CMAME190

            C1 = mateprop(6);
            mu = 2*C1;
            
        dpmat3 = [8 0 0 
                  0 0 0 
                  0 0 2 
                  0 0 0 
                  0 8 0 
                  0 0 2 
                  0 0 2 
                  0 0 2 
                  2 2 0 ];

            dmati=mu*dpmat3;
            
        end
    end
end