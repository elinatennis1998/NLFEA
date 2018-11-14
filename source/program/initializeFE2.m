% Tim Truster
% 04/19/2013
%
% Subroutine to initialize FE arrays for NLFEA ver2 in accordance with the
% FEAP standard

FE_err = 0;

%     Allocate and zero arrays

idFEAP2 = zeros(ndf,numnp);
idFEAPBC2 = idFEAP2;
idFEAPBCnp2 = idFEAP2;

%     Input boundary conditions

pbouin2

%       Perform check on mesh again to set final boundary codes

pidset2

%     Determine current profile, set current profile

seteq2

%     Set flags
iter = 0;
step = 0;

%     Input nodal and surface loads

ploadi2

%     Set body force internal loads

if(abs(intfl2))
pbodyf2
end

% % Read and interpret boundary conditions, loads, (initial conditions)
% 
% assign_bc_load_dataNL
