% Tim Truster
% 01/04/2012
% UIUC

% Master switch routine for setting which element subroutine is called
% Called by: FormFE, assign_bc_load_dataNL, FormS

if exist('nel','var') == 0
    nel = 4;
end

% Zero-out the proper elemental assembly arrays
switch isw %Task Switch
    
    case -1 % Surface Loads
        ElemF = zeros(nst,1);
    case 1 %Get Material Properties
        
    case 3 %Get Stiffness, Force
        ElemF = zeros(nst,1);
        ElemK = zeros(nst,nst);
        ElemFn = zeros(nst,1);
    case 5 %Get Mass
        ElemM = zeros(nst,nst);
    case 6 %Get Force
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);
    case 9 %Get global error
        ElemF = zeros(nst,1);
    case 11 %Get Error Indicators
        ElemE = zeros(numEn,1);
    case 12 % system energy
        ElemE = 0;
    case 15 %Body forces
        ElemF = zeros(nst,1);
    case 16 % J-Integral
        ElemJ = zeros(ndm,1);
    case 21 %Get Stiffness
        ElemK = zeros(nst,nst);
    case 22 % L2 stress projection
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);
    case 25 % Stress averaging
        ElemS = zeros(nel,npstr+1);
    case 26 % Elemental stresses
        ElemS = zeros(nel,nestr);
    case 51% Volume stress/strain homogenization
        ElemSS = zeros(13,1);
    case 98% Lie group calculation, end of step
        ElemL = zeros(nel,numLie+1);
    case 100% Lie group calculation, during step
        ElemL = zeros(nel,numLie+1);
end

	if nonlin == 0 %linear analysis
        switch iel
            case 1 %Small-Deformation Isotropic Elastostatics Element
                if ndm == 3
                    NL_Elem8_3d % L_Elem1_3d % 
                else %ndm == 2
                    L_Elem1_2d % L_Elem1_2d2 % 
                end
            case 2 % Poisson Equation
                if ndm == 3

                else %ndm == 2
                    L_Elem2_2dVMS
                end
            case 3 %Stabilized Mixed Pressure-Displacement Element
                if ndm == 2
                    L_Elem3_2dVMS % L_Elem3_2dVMSP1P0 % 
%                     NL_Elem3_2d
                else %ndm == 3
                    L_Elem3_3d2
                end
            case 4 %Small-Deformation Mixed Press-Displac DG element
                if ndm == 2
                    L_Elem3_2dDG % L_Elem3_2dDGP1P0 % 
                else %ndm == 3

                end
            case 5 %Small-Deformation DG element, with triangular segments
                if ndm == 2
                    L_Elem1_2dDG
                else %ndm == 3
                    L_Elem1_3dDG2 % L_Elem1_3dDG % 
                end
            case 6 % Mixed DG element with u-p one side, u-only on other
                if ndm == 2
%                     NL_Elem5_2d2
                    L_Elem6_2dDG
                else %ndm == 3

                end
            case 7 % Darcy Flow interior
                if ndm == 2
                    L_Elem4_2dVMSb
                else %ndm == 3

                end
            case 8 %Small-Deformation DG element, full element sectors
                if ndm == 2
                    L_Elem8_2dDG
                else %ndm == 3
                    L_Elem8_3dDG
                end
            case 9 % 3d DG element in normal direction, linear spring in tangential
                if ndm == 3
                    NL_Elem8_3dDG % NL_Elem8_3dDG2 % 
                else %ndm == 2
                    NL_Elem9_2d % multiscale dynamics; does not work
                end
            case 10 % Darcy Flow DG
                if ndm == 2
                    L_Elem4_2dDGb
                else %ndm == 3

                end
            case 11 % Darcy-Stokes
                if ndm == 3
                    
                else %ndm == 2
                    NL_Elem11_2dVMS
                end
            case 12 % Darcy Flow DG
                if ndm == 2
                    L_Elem6_2dDG3b
                else %ndm == 3

                end
            case 13 % Darcy-Stokes Flow DG
                if ndm == 2
                    L_Elem13_2d
                else %ndm == 3

                end
            case 14 % Axial rod
                if ndm == 1
                    L_Elem0_1d
                elseif ndm == 2
                    L_Elem0_2d
                else %ndm == 3

                end
            case 15 % 1d beam
                if ndm == 1
                    L_Elem7_1d
                elseif ndm == 2
                    L_Elem7_2d
                else %ndm == 3

                end
            case 16 %Stabilized Mixed Pressure-Displacement Element
                if ndm == 2
                    L_Elem3_2dVMSP1P0 % L_Elem3_2dVMS % 
%                     NL_Elem3_2d
                else %ndm == 3
                    L_Elem3_3d2
                end
            case 17 %Implicit Error Element
                if ndm == 2
                    L_Elem3_2dDGP1P0 % L_Elem3_2dDG % 
                else %ndm == 3

                end
            case 18 %Implicit Error Element
                if ndm == 2
                    L_Elem3_2dDGP1P0_nonconform % L_Elem3_2dDG % 
                else %ndm == 3

                end
            case 19 %Vectorized Poisson Equation - Acoustic Tensor for Finite Deformations
                if ndm == 2
                    L_Elem19_2d
                else %ndm == 3

                end
            case 21 % Generalized alpha method dynamics
                if ndm == 2
                    NL_Elem5_2d2
                else %ndm == 3
                    L_Elem5_3d
                end
            case 22 % bi-linear cohesive zone element; element #7 in DEIP
                if ndm == 2
                    L_Elem22_2dDG
                else % ndm == 3
                    L_Elem22_3dDG
                end
            case 23 %Small-Deformation periodic boundary conditions using Lagrange multipliers
                L_Elem23_3d
            case 24 %Small-Deformation multi-point constraint using Lagrange multipliers
                L_Elem24_3d
            case 25 %Penalty element
                L_Elem25_3d
            case 26 % Three-field Domain Decomposition element
                if ndm == 3
                    if nen == 20 % only 2 constituents, 8 node hexahedral for solid elements
                    NL_Elem34_3d
                    else % 3 constituents
                    NL_Elem34_3dmulti
                    end
                else %ndm == 2

                end
            case 27 % Reaction-Diffusion element
                if ndm == 1
                    L_Elem27_1dRFB
                else %ndm == 3

                end
            case 28 % Reaction-Diffusion element with tau
                if ndm == 1
                    L_Elem28_1d
                else %ndm == 3

                end
            case 31 % Nitsche plane-stress 2 bar coupling
                if ndm == 2
                    L_Elem30_1d
                else %ndm == 3

                end
            case 32 % Nitsche plane-stress 2 beam coupling
                if ndm == 2
                    L_Elem32_1d
                else %ndm == 3

                end
            case 36 % Helmholtz element
                if ndm == 1
                    L_Elem36_1d
                else %ndm == 2
                end
            case 37 % 1d wave element for GLS-GGLS method %Harari's actual method
                if ndm == 1
                    L_Elem37_1dGGLS % L_Elem37_1dHarari % 
                else %ndm == 2
                end
            case 38 % 1d wave equation average mass, VMS
                if ndm == 1
                    NL_Elem5_1dEquiv % NL_Elem5_1d % 
                else %ndm == 2
                end
            case 39 % 1d dynamic rod element VMS, expanded matrix version
                if ndm == 1
                    NL_Elem39_1d2 % NL_Elem39_1d % 
                else %ndm == 3

                end
            case 51 % 2D-Truss (Omar Nassif)
                L_Omar_2dtruss
            case 52 % 3D-Truss (Omar Nassif)
                L_Omar_3dtruss
            case 53 % 2D-Frame (Omar Nassif)
                L_Omar_2dframe
            case 56 %2D-Bar (Elina Geut)
                L_Elem56_2d_Elina
            case 57 %1-D DG linear element (Elina Geut)
                L_Elem57_1dDG
            case 59 %1D bar linear element (Elina Geut)
                if exist('pL','var') && pR ~= pL
                    L_Elem59_1d_2_Elina
                else
                    L_Elem59_Elina
                end
            case 60 %1D Quadratic
                L_Elem60_Elina
            case 61 % 1D bar linear element
                if ndm ==1
                    L_Elem61_1d
                end 
            case 62 % 1D bar linear DG element
                if ndm == 1
                    L_Elem62_1dDG
                end
            otherwise
                error('Element type not implemented')
        end
    else %nonlinear analysis
        switch iel
            case 1 %Nodal smooth/rough contact
                if ndm == 3
                    NL_Elem1_3d
                else %ndm == 2
                    NL_Elem1_2d
                end
            case 2 % Pure displacement
                if ndm == 3
                    NL_Elem2_3dM % NL_Elem2_3d % 
                else %ndm == 2
                    if exist('use_function','var') == 1 && use_function == 1
                    NL_Elem2_2dF % function-based version
                    else
                    NL_Elem2_2dM % script-based version
                    end
                end
            case 3 %Linear Elastodynamics
                if ndm == 3
                else %ndm == 2
                    NL_Elem3_2d
                end
            case 4 %Unequal tension/shear
                if ndm == 3
                    NL_Elem4_3d
                else %ndm == 2
                    NL_Elem4_2d
                end
            case 5 % Mixed form
                if ndm == 3
                    if nel == 4
                        NL_Elem5_3dNSCST % NL_Elem5_3dSCST % NL_Elem5_3dCST % NL_Elem5_3dCST2 % 
                    else
                        NL_Elem5_3dNS % NL_Elem5_3dS % NL_Elem5_3dM % 
                    end
                else %ndm == 2
                    if nel == 3
                        NL_Elem5_2dNSCST % NL_Elem5_2dNS % NL_Elem5_2dSCST % NL_Elem5_2dS % NL_Elem5_2dCST % 
                    else
                        NL_Elem5_2dM  % NL_Elem5_2dNS % NL_Elem5_2d4 % NL_Elem5_2dS %  NL_Elem5_2dFS % NL_Elem5_2d2 % NL_Elem5_2d2ns % NL_Elem5_2d2h % 
                    end
                end
            case 6 % Generalized Alpha method
                if ndm == 3
                    NL_Elem6_3d3
                    % Use NL_Elem6_3d3 for simulation, NL_Elem6_3d for FormI
                else %ndm == 2
                    NL_Elem5_2d2
                end
            case 7 % Energy constraint small defo
                if ndm == 3
                    NL_Elem7_3d
                else %ndm == 2
                    NL_Elem7_2d
                end
            case 8 % Energy-momentum conserving
                if ndm == 3
                    NL_Elem8_3d
                else %ndm == 2
                    NL_Elem8_2d
                end
            case 9 %Mixed Integral Coulomb
                if ndm == 3
                    NL_Elem8_3dDG
                else %ndm == 2
                    NL_Elem9_2d4 % NL_Elem9_2d3 % NL_Elem9_2d2 % NL_Elem9_2d % 
                end
            case 10 % Darcy Flow DG
                if ndm == 2
                    NL_Elem10_2d2 % NL_Elem10_2d % 
                else %ndm == 3

                end
            case 11 %Constitutive Model
                if ndm == 3

                else %ndm == 2
                    NL_Elem10_2d
                end
            case 12 %Constitutive Model
                if ndm == 3

                else %ndm == 2
                    NL_Elem12_2d
                end
            case 13 % Bbar small strain J2 plasticity
                if ndm == 3
                    if nen == 10
                    Bbar3d_Elem4 % 
                    else
                    Bbar3d_Elem % Bbar3d_Elem2 % 
                    end
                else %ndm == 2
                    if nen == 6
                    Bbar2d_Elem4 % 
                    else
                    Bbar2d_Elem3 % Bbar2d_Elem2 % Bbar2d_Elem % 
                    end
                end
            case 14
                if ndm == 3
                    NL_Elem2_3dFbar % NL_Elem2_3dMIE
                else %ndm == 2
                    NL_Elem2_2dMCPFbar % NL_Elem2_2dMCP % NL_Elem2_2dFbar % NL_Elem2_2dMIE % 
                end
            case 15
                if ndm == 3
                else %ndm == 2
                    NL_Elem11_2d
                end
            case 16
                if ndm == 3
                    
                else %ndm == 2
                    NL_Elem11_2dDG
                end
            case 17
                if ndm == 3
                    
                else %ndm == 2
                    NL_Elem2_2dMDG
                end
            case 18
                if ndm == 3
                    if igrow == 1
                    NL_Elem2_3dMG % NL_Elem2_3dMGc2 % NL_Elem2_3dMGc % 
                    else
                    NL_Elem2_3dMGc2 % NL_Elem2_3dMG % NL_Elem2_3dMGc % 
                    end
                else %ndm == 2
                end
            case 19
                if ndm == 3
                    NL_Elem2_3dMGd
                else %ndm == 2
                end
            case 20
                if ndm == 2
                    NL_Elem20_2d
                else %ndm == 2
                end
            case 21  % large deformation 
                if ndm == 2
                    NL_Elem21_2d_2 %NL_Elem21_2d_4 %NL_Elem21_2d_5 %   %    NL_Elem21_2d_3 %  %  
                else %ndm == 2
                    NL_Elem21_3d_7 % NL_Elem21_3d_8 % NL_Elem2_3dMDG % 
                end 
            case 22   %debond case
                if ndm == 2
                 
                else %ndm == 2
                   NL_Elem21_3d_8 % 
                end                 
            case 23 %weakly enforced boundary
                if ndm == 2
                    NL_Elem23_2d_1
                else %ndm == 2                   
                end     
            case 24 %mixed method
                if ndm == 2
                    NL_Elem24_2d_1
                else %ndm == 2                   
                end                 
            case 25
                if ndm == 2
                    NL_Elem25_2d
                else %ndm == 2
                end
            case 26 % Density growth, mixture theory
                if ndm == 3
                    NL_Elem26_3dMG
                else %ndm == 2
                end
            case 27 %31
                if ndm == 3
                    NL_Elem2_3dMGd_2con
                else %ndm == 2
                end
            case 28 %32
                if ndm == 3
                    NL_Elem2_3dMGc_2con
                else %ndm == 2
                end
            case 29 %33
                if ndm == 3
                    NL_Elem33_3dM
                else %ndm == 2
                end
            case 30 %34
                if ndm == 3
                    NL_Elem34_3dBoth
                else %ndm == 2
                end
            case 31 % Nitsche plane-stress 2 bar coupling
                if ndm == 2
                    L_Elem30_1d
                else %ndm == 3

                end
            case 32 % Nitsche plane-stress 2 beam coupling
                if ndm == 2
                    L_Elem32_1d
                else %ndm == 3

                end
            case 34 % 1d dynamic rod element
                if ndm == 1
                    NL_Elem0_1dGLS % NL_Elem0_1d % 
                else %ndm == 2
                end
            case 35 % 1d dynamic rod element VMS
                if ndm == 1
                    NL_Elem2_1dRFB % NL_Elem2_1d % 
                else %ndm == 2
                end 
            case 36 % 1d wave equation Harari
                if ndm == 1
                    NL_Elem5_1d2 % NL_Elem5_1dEquiv % NL_Elem5_1d % 
                else %ndm == 2
                end
            case 37 % 1d wave equation RFB
                if ndm == 1
                    L_Elem37_1dGGLS % L_Elem37_1dHarari % NL_Elem5_1dRFB2 % NL_Elem5_1dRFB % 
                else %ndm == 2
                end
            case 38 % 1d wave equation average mass, VMS
                if ndm == 1
                    NL_Elem5_1dEquiv % NL_Elem5_1d % 
                else %ndm == 2
                end
            case 39 % 1d wave equation average mass, VMS with history as nodes
                if ndm == 1
                    NL_Elem39_1d2 % NL_Elem39_1d % 
                else %ndm == 2
                end
            case 40 % Generalized alpha method dynamics - linear kinematics element, but nonlinear algorithmic implementation
                if ndm == 2
                    NL_Elem5_2dNL
                else %ndm == 3

                end
            case 41 % qausi-static composite debonding
                if ndm == 2
                    NL_Elem41_2dDG2 % NL_Elem41_2dDG % 
                else %ndm == 3
                    NL_Elem3_3d_bubble % NL_Elem3_3d3chatter % NL_Elem4_3d % NL_Elem3_3d3triang % 
                end
            case 42 % bi-linear cohesive zone element
                if ndm == 2
                    NL_Elem42_2dDG2
                else
                    L_Elem7_3dDG
                end
            case 46 % Hypoelasto-plasticity from WARP3D
                if ndm == 2
                    NL_Elem46_2d
                else %ndm == 3
                    NL_Elem46_3d
                end
            case 47 % Mark CP version 1 from WARP3D
                if ndm == 2
                    NL_Elem47_2d
                else %ndm == 3
                    NL_Elem47_3d
                end
            case 48 % Explicit integrator for CP, 1-pt rule
                if ndm == 2
                    error('')
                else %ndm == 3
                    NL_Elem48_3d
                end
            case 51 % 2d von Mises plasticity pure-displacement
                if ndm == 2
                    NL_Elem51_2d
                end
            case 52 % 2d von Mises plasticity Q2Q1 u-p mixed formulation
                if ndm == 2
                    NL_Elem52_2d
                elseif ndm == 3
                    NL_Elem52_3d
                end
            case 53 % 2d von Mises plasticity pure-displacement DG
                if ndm == 2
                    NL_Elem53_2d
                end
            case 54 % 2d von Mises plasticity Q2Q1 u-p mixed DG
                if ndm == 2
                    NL_Elem54_2d
                elseif ndm == 3
                    NL_Elem54_3d
                end
            case 55 % 2d von Mises plasticity coupled DG
                if ndm == 2
                    NL_Elem55_2d
                end
            case 56 % 1d FDM element, starting as ADR element
                if ndm == 1
                    NL_Elem56_1d
                end
            case 57 % template 1d element with 3-component solution field
                if ndm == 1
                    NL_Elem57_1dVMS % NL_Elem57_1dLS % NL_Elem57_1d % 
                end
            case 61 % 1D bar nonlinear Element
                if ndm == 1
                    NL_Elem61_1d
                end
               
            case 62 % 1D bar nonlinear DG element
                if ndm == 1
                   NL_Elem62_1dDG
                end 
            case 63 % 1D bar nonlinear DG element(general)
                if ndm == 1
                   NL_Elem63_1dDGF
                   %NL_Elem63_1dDGzeta100
                   %NL_Elem63_1dDGzetachange
                end 
            otherwise
                error('Element type not implemented')
        end
	end
