% Tim Truster
% 04/19/2013
%
% Subroutine to initialize FE arrays for NLFEA ver2 in accordance with the
% FEAP standard

FE_err = 0;

%     Set pointers for allocation of mesh arrays

nie       = ndf + 8;

%     Allocate and zero arrays

ieFEAP = zeros(nie,nummat);
iedof = zeros(ndf,nen,nummat);
ixFEAP = zeros(7,numel);
idFEAP = zeros(ndf,numnp);
idFEAPBC = idFEAP;
idFEAPBCnp = idFEAP;

%     Input boundary conditions

pbouin

%     Data input for material sets

pmatin

%       Perform check on mesh again to set final boundary codes

pidset

%     Set up stress history addresses

sethis

%     Determine current profile, set current profile

seteq

nha = 1;
nhb = nhmax + 1;
nhc = 2*nhmax + 1;

%     Set values for nels, other flags

psflags

%     Input nodal and surface loads

ploadi

%     Set body force internal loads

if(abs(intfl))
pbodyf
end
