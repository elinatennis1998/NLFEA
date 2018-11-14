% Tim Truster
% 04/19/2013
%
% Subroutine to initialize FE arrays for NLFEA ver2 in accordance with the
% FEAP standard

FE_err = 0;

%     Allocate and zero arrays

idFEAP1 = zeros(ndf,numnp);
idFEAPBC1 = idFEAP1;
idFEAPBCnp1 = idFEAP1;

%     Input boundary conditions

pbouin1

%       Perform check on mesh again to set final boundary codes

pidset1

%     Determine current profile, set current profile

seteq1

%     Set flags
iter = 0;
step = 0;

%     Input nodal and surface loads

ploadi1

%     Set body force internal loads

if(abs(intfl1))
pbodyf1
end

% % Read and interpret boundary conditions, loads, (initial conditions)
% 
% assign_bc_load_dataNL
