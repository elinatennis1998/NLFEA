
% Tim Truster
% 12/09/2013
% 
% Program to generate nodal loads for finite element analysis.
% Input file is selected interactively.
%
% Nodal loads that form the external force vector are calculated using
% exactly the same functions that are used in NL_FEA_Program; this allows
% the user to verify that the input file and element routine are
% functioning as intended.
%
% Loads are stored in FcMat and FcMatnp for proportional and
% nonproportional loads, respectively.
%
% Boundary conditions are ignored so that all loads can be computed.
%
% To ignore body forces, remove those flags from the input file.
%
% These values are suitable for writing to a FEAP input file.

%% Ask for input file name and load the file
if exist('filename.mat','file')
clear load
load('filename.mat', 'defaultname','defaultpath')
else
defaultname = 'DNE';
end
fprintf('Provide name of input file (default is %s)\n',defaultname);
[filename,pathname] = uigetfile('*.m','Select the NLFEA input file');
if isnumeric(filename)
     pathname = defaultpath;
     filename = defaultname;
end
defaultname = filename;
defaultpath = pathname;
save('filename.mat','defaultname','defaultpath');
run([pathname filename(1:end-2)])

% Read major problem size parameters and solution options
pstart

% Initialize default parameters
pdefault

% Set pointers for allocation of mesh arrays
nie = ndf + 8;

% Allocate and zero arrays
ieFEAP = zeros(nie,nummat);
iedof = zeros(ndf,nen,nummat);
ixFEAP = zeros(7,numel);
ixFEAP(7,:) = RegionOnElement(:);

% Allocate full-size NDOFT without boundary conditions
neq = numnp*ndf;
NDOFT = reshape((1:numnp*ndf),numnp,ndf);

% Data input for material sets
pmatin

% Set values for nels, other flags
psflags

% Input nodal and surface loads
ploadi

% Set body force internal loads
if(abs(intfl))
pbodyf
end

% Reshape the load vector into a matrix
FcMat = reshape(Fc1,numnp,ndf);
FcMatnp = reshape(Fc1np,numnp,ndf);
FcMatU = reshape(FcU,numnp,ndf);