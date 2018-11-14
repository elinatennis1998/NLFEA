% Default parameter values
% Tim Truster
% 10/06/2013

%% Default values for solution algorithm parameters

% Load algorithm library data
palgolib
    
switch transient
    
    case -1 % Linear static analysis
        if ~exist('s_del_a','var')
        s_del_a = 1;
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
    case 0 % Quasi-static Proportional loading
        
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('mults','var')
        mults = 1; % values of proportional load factor
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('smax','var')
        smax = 1; % maximum scaling parameter for line-search
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('Residratio','var')
        Residratio = 10^-11; % relative tolerance for residual norm
        end
%         if ~exist('ddmax','var')
%         ddmax = 0; % ?
%         end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        if ~exist('lsOn','var') % Line search flag; true for on
        lsOn = 0;
        end
    case {1,2} % Dynamic HHT-Alpha Method
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
        % Currently only works for NL_Elem3_2d and NL_Elem7_2d
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for energy conserving
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 4
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 5 % Generalized alpha method
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
%         rhoinf = 0.75;
%         Nalphaf = rhoinf/(rhoinf  + 1);0;(1 - rhoinf)/(rhoinf + 1);
%         Nalpham = (2*rhoinf - 1)/(rhoinf + 1);0;(rhoinf - 1)/(1 + rhoinf);
%         Nbeta = 1/4*(1 - Nalpham + Nalphaf)^2;1/4;
%         Ngamma = 1/2 - Nalpham + Nalphaf;1/2;
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('Nalphaf','var')
        Nalphaf = 0; % time-integration parameter
        end
        if ~exist('Nalpham','var')
        Nalpham = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 6 % Multiscale Dynamics
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('Nalphaf','var')
        Nalphaf = 0; % time-integration parameter
        end
        if ~exist('Nalpham','var')
        Nalpham = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 7 % Quasi-static Arc-Length Method
        if ~exist('b','var')
        b = 1; % [0 1] = [disp control, force control]
        end
        if ~exist('s_del_a','var')
        s_del_a = 1; % measure of arc-length distance
        end
        if ~exist('lamdamax','var')
        lamdamax = 1; % maximum load-level parameter
        end
        if ~exist('smax','var')
        smax = 1; % maximum scaling parameter for line-search
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-11; % relative tolerance for residual norm
        end
        if ~exist('residratio','var')
        residratio = 10^-12; % relative tolerance for arc-length norm
        end
        if ~exist('ddmax','var')
        ddmax = 0; % ?
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equalting multiples of this #
        end
        
    case 8 % Thermoelasticity with Newmark for disp and gener-trap for temp
        if ~exist('enercons','var')
        enercons = 0; % flag=1 for printing energy norm
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('Nalphaf','var')
        Nalphaf = 0; % time-integration parameter for displacement
        end
        if ~exist('Nalpham','var')
        Nalpham = 0; % time-integration parameter for displacement
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter for displacement
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter for displacement
        end
        if ~exist('Nalphat','var')
        Nalphat = 1; % time-integration parameter for temperature
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equalting multiples of this #
        end
        
    case 10 % Explicit Dynamics
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 0; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 3; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 11 % Crystal plasticity with averaging of Fe_inv
        
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('mults','var')
        mults = 1; % values of proportional load factor
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('smax','var')
        smax = 1; % maximum scaling parameter for line-search
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('Residratio','var')
        Residratio = 10^-11; % relative tolerance for residual norm
        end
%         if ~exist('ddmax','var')
%         ddmax = 0; % ?
%         end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        if ~exist('numsubcycles','var')
        numsubcycles = 10; % dump data at steps equaling multiples of this #
        end
        if ~exist('subcycTol','var')
        subcycTol = 1e-8;1e-3; % dump data at steps equaling multiples of this #
        end
        
    case 12 % Forward transient Euler, based off of case 10
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 0; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 3; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 13 % GLS
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for energy conserving
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('Nbeta','var')
        Nbeta = 1/4; % time-integration parameter
        end
        if ~exist('Ngamma','var')
        Ngamma = 1/2; % time-integration parameter
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
    case 14 % generalized trapezoidal algorithm, based off of case 5
        if ~exist('enercons','var')
        enercons = 1; % flag=1 for printing energy norm
        end
        if ~exist('tstep','var')
        tstep = 1; % constant time step size
        end
        if ~exist('Nalpha','var')
        Nalpha = 0; % time-integration parameter
        end
        if ~exist('itermax','var')
        itermax = 10; % maximum number of iterations
        end
        if ~exist('stepmax','var')
        stepmax = 1; % number of steps
        end
        if ~exist('datastep','var')
        datastep = stepmax; % size of output arrays
        end
        if ~exist('Residratio','var')
        Residratio = 10^-12; % relative tolerance for residual norm
        end
        if ~exist('reststep','var')
        reststep = 100; % dump data at steps equaling multiples of this #
        end
        
end

% BC/IC/Load input counters, NCR=1, and error checks
if ~exist('numComp','var')
    numComp = 0; % Tying node commands
elseif numComp > 0 && ~exist('NodeComp','var')
    error('NodeComp is not defined')
elseif numComp > 0 && numComp < size(NodeComp,1)
    disp('Warning, node tying: numComp < size(NodeComp,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numBC','var')
    numBC = 0; % Boundary condition, scaled by lamda
elseif numBC > 0 && ~exist('NodeBC','var')
    error('NodeBC is not defined')
elseif numBC > 0 && numBC < size(NodeBC,1)
    disp('Warning, node BC: numBC < size(NodeBC,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numBCnp','var')
    numBCnp = 0; % Boundary condition, non-proportional
elseif numBCnp > 0 && ~exist('NodeBCnp','var')
    error('NodeBCnp is not defined')
elseif numBCnp > 0 && numBCnp < size(NodeBCnp,1)
    disp('Warning, node BCnp: numBCnp < size(NodeBCnp,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numBCmt','var') % tabular BCs
    numBCmt = 0;
elseif numBCmt > 0 && ~exist('NodeBCmt','var')
    error('NodeBCmt is not defined')
end

if ~exist('numIC','var')
    numIC = 0; % Initial field value (e.g. displacements)
elseif numIC > 0 && ~exist('NodeIC','var')
    error('NodeIC is not defined')
elseif numIC > 0 && numIC < size(NodeIC,1)
    disp('Warning, node uIC: numIC < size(NodeIC,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numVIC','var')
    numVIC = 0; % Initial rates (e.g. velocities)
elseif numVIC > 0 && ~exist('NodeVIC','var')
    error('NodeVIC is not defined')
elseif numVIC > 0 &&numVIC < size(NodeVIC,1)
    disp('Warning, node vIC: numBC < size(NodeVIC,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numNodalF','var')
    numNodalF = 0; % Nodal loads, scaled by lamda
elseif numNodalF > 0 && ~exist('NodeLoad','var')
    error('NodeLoad is not defined')
elseif numNodalF > 0 && numNodalF < size(NodeLoad,1)
    disp('Warning, node loads: numNodalF < size(NodeLoad,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numNodalFnp','var')
    numNodalFnp = 0; % Nodal loads, non-proportional
elseif numNodalFnp > 0 && ~exist('NodeLoadnp','var')
    error('NodeLoadnp is not defined')
elseif numNodalFnp > 0 && numNodalFnp < size(NodeLoadnp,1)
    disp('Warning, node loads: numNodalFnp < size(NodeLoadnp,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numSL','var')
    numSL = 0; % Surface loads, scaled by lamda, isw=-1
elseif numSL > 0 && ~exist('SurfacesL','var')
    error('SurfacesL is not defined')
elseif numSL > 0 && numSL < size(SurfacesL,1)
    disp('Warning, surface loads: numSL < size(SurfacesL,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numSLnp','var')
    numSLnp = 0; % Surface loads, non-proportional, isw=-1
elseif numSLnp > 0 && ~exist('SurfacesLnp','var')
    error('SurfacesLnp is not defined')
elseif numSLnp > 0 && numSLnp < size(SurfacesLnp,1)
    disp('Warning, surface loads_np: numSLnp < size(SurfacesLnp,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numSLmt','var')
    numSLmt = 0; % Surface loads, scaled by lamda, isw=-1
elseif numSLmt > 0 && ~exist('SurfacesLmt','var')
    error('SurfacesLmt is not defined')
end

if ~exist('numD','var')
    numD = 0; % Weak Dirichlet BCs, currently inactive, isw=-3
end

if ~exist('intfl','var')
    intfl = 0; % T/F flag to include body forces
elseif ~exist('numBF','var')
    numBF = 0; % Number of material sets with body forces
elseif numBF > 0 && ~exist('BodyForce','var')
    error('BodyForce is not defined')
elseif numBF > 0 && numBF < size(BodyForce,1)
    disp('Warning, body forces: numBF < size(BodyForce,1)')
    disp('Paused: press any key to continue')
end


if ~exist('Uiterstep','var')
    Uiterstep = 0; % Evaluate user loads at 0=step, 1=every iteration
end

if ~exist('numUSL','var')
    numUSL = 0; % Surface loads defined by user at each step/iteration
elseif numUSL > 0 && ~exist('SurfacesU','var')
    error('SurfacesU is not defined')
elseif numUSL > 0 && numUSL < size(SurfacesU,1)
    disp('Warning, surface loadsU: numUSL < size(SurfacesU,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numUBF','var')
    numUBF = 0; % Body force defined by user at each step/iteration
elseif numUBF > 0 && ~exist('BodyForceU','var')
    error('BodyForceU is not defined')
elseif numUBF > 0 && numUBF < size(BodyForceU,1)
    disp('Warning, body forcesU: numUBF < size(BodyForceU,1)')
    disp('Paused: press any key to continue')
end


% Load input counters, NCR=[2,3]
if ~exist('numSL_new','var')
    numSL_new = 0; % Surface loads, scaled by lamda, isw=-1
elseif numSL_new > 0 && ~exist('SurfacesL_new','var')
    error('SurfacesL_new is not defined')
elseif numSL_new > 0 && numSL_new < size(SurfacesL_new,1)
    disp('Warning, new surface loads: numSL_new < size(SurfacesL_new,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numSLnp_new','var')
    numSLnp_new = 0; % Surface loads, non-proportional, isw=-1
elseif numSLnp_new > 0 && ~exist('SurfacesLnp_new','var')
    error('SurfacesLnp_new is not defined')
elseif numSLnp_new > 0 && numSLnp_new < size(SurfacesLnp_new,1)
    disp('Warning, new surface loads_np: numSLnp_new < size(SurfacesLnp_new,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numNodalF_new','var')
    numNodalF_new = 0; % Nodal loads, scaled by lamda
elseif numNodalF_new > 0 && ~exist('NodeLoad_new','var')
    error('NodeLoad_new is not defined')
elseif numNodalF_new > 0 && numNodalF_new < size(NodeLoad_new,1)
    disp('Warning, new nodal loads: numNodalF_new < size(NodeLoad_new,1)')
    disp('Paused: press any key to continue')
end

if ~exist('numNodalFnp_new','var')
    numNodalFnp_new = 0; % Nodal loads, non-proportional
elseif numNodalFnp_new > 0 && ~exist('NodeLoadnp_new','var')
    error('NodeLoadnp_new is not defined')
elseif numNodalFnp_new > 0 && numNodalFnp_new < size(NodeLoadnp_new,1)
    disp('Warning, new nodal loads_np: numNodalFnp_new < size(NodeLoadnp_new,1)')
    disp('Paused: press any key to continue')
end


if ~exist('numSI','var')
    numSI = 0; % Interface elements
elseif numSI > 0 && ~exist('SurfacesI','var')
    disp('SurfacesI is not defined; might be okay')
    disp('Paused: press any key to continue')
elseif numSI > 0 && numSI < size(SurfacesI,1)
    disp('Warning, interfaces: numSI < size(SurfacesI,1)')
    disp('Paused: press any key to continue')
end

% Error estimation macro
if ~exist('expliciterr','var')
expliciterr = 0; % T/T flag for performing explicit error estimation
end
if ~exist('impliciterr','var')
impliciterr = 0; % T/T flag for performing implicit and global error estimation
end
if ~exist('implicon','var')
implicon = 1; % T/T flag for using local-implicit or local-explicit error for driving global error
end
if ~exist('subm','var')
subm = 1; % T/T flag for ?
end
if ~exist('residflag','var')
residflag = 0; % T/T flag for ?
end

% Dynamics flags
if ~exist('lumping','var')
lumping = 0; % set 0 for no mass lumping, 1 for lumping of mass matrix
end

% Algorithm parameters
if ~exist('restartmat','var')
restartmat = 1; % T/T flag for writing restart (.mat) files
end
if ~exist('printRnorm','var')
printRnorm = 1; % T/T flag for printing residual norm to screen
end
if ~exist('hrstore','var')
hrstore = 0; % T/T flag for storing history variables for all time levels
end

% Other parameters
if ~exist('nelPn','var')
% Mixed element and stress nodes per element
nelPn = 1; %1 for unequal-order, 0 for equal-order
end
if ~exist('LieLump','var')
LieLump = 1; % flag for doing lumped mass projection or just extrapolation for Lie group
end
