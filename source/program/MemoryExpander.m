% Tim Truster
% 09/26/2014
%
% Memory expander for NL_FEA_Program output files
% Run this when restarting or continuing a simulation, etc to add more data
% rows to the history arrays.
%
% Note: additional arrays can be defined for allocating stresses, etc for 
% not every step.
%
% NOTE: Currently, for restar analyses, the restart '.mat' file should bec
% modified first using this routine and then re-saved, because NCR=3 in
% NL_FEA_Program clears the workspace before loading the restart data, so
% this file fails at that point. Instead, increase the sizes and save a
% file outside of running NL_FEA.

% Set a default value for number of additional steps
if ~exist('datastep_add','var') 
    datastep_add = 10;
end
datastep_new = datastep + datastep_add;

% Displacements
if exist('DispList','var')
    DispListOld = DispList;
    DispList = zeros(ndf,numnp,datastep_new);
    DispList(:,:,1:datastep) = DispListOld;
    clear('DispListOld')
end

% Velocities
if exist('VeloList','var')
    VeloListOld = VeloList;
    VeloList = zeros(ndf,numnp,datastep_new);
    VeloList(:,:,1:datastep) = VeloListOld;
    clear('VeloListOld')
end

% Accelerations
if exist('AcceList','var')
    AcceListOld = AcceList;
    AcceList = zeros(ndf,numnp,datastep_new);
    AcceList(:,:,1:datastep) = AcceListOld;
    clear('AcceListOld')
end

% Reactions
if exist('ForcList','var')
    ForcListOld = ForcList;
    ForcList = zeros(ndf,numnp,datastep_new);
    ForcList(:,:,1:datastep) = ForcListOld;
    clear('ForcListOld')
end

% Crystal plasticity
if exist('LieHist','var') && LieHist == 1
    LieHListOld = LieHList;
    LieHList = zeros(numLie,numnp,datastep_new);
    LieHList(:,:,1:datastep) = LieHListOld;
    clear('LieHListOld')
    SubnormsTOld = SubnormsT;
    SubnormsT = zeros(numSubnormvals,itermax+1,datastep_new);
    SubnormsT(:,:,1:datastep) = SubnormsTOld;
    clear('SubnormsTOld')
    SubnormsLOld = SubnormsL;
    SubnormsL = zeros(numSubnormvals,(itermax+1)*(numsubcycles+1),datastep_new);
    lenStre = size(SubnormsLOld,2);
    SubnormsL(:,1:lenStre,1:datastep) = SubnormsLOld(:,1:lenStre,1:datastep);
    clear('SubnormsLOld')
end

% Nodal Stresses
if SHist
    if ~exist('SHListAdd','var')
        SHListAdd = 1:datastep_add;
        SHListAdd = stepmax + SHListAdd;
    end
    SHList = [SHList SHListAdd];
    lenStre = size(StreList,3);
    StreListOld = StreList;
    StreList = zeros(npstr,numnp,lenStre+datastep_add);
    StreList(:,:,1:lenStre) = StreListOld;
    clear('StreListOld')
end

% Element Stresses
if SEHist
    if ~exist('SEHListAdd','var')
        SEHListAdd = 1:datastep_add;
        SEHListAdd = stepmax + SEHListAdd;
    end
    SEHList = [SEHList SEHListAdd];
    lenStre = size(StreListE,3);
    StreListEOld = StreListE;
    StreListE = zeros(nestr,numel,lenStre+datastep_add);
    StreListE(:,:,1:lenStre) = StreListEOld;
    clear('StreListEOld')
end

% Explicit error
if expliciterr
    if ~exist('EEHListAdd','var')
        EEHListAdd = 1:datastep_add;
        EEHListAdd = stepmax + EEHListAdd;
    end
    EEHList = [EEHList EEHListAdd];
    lenStre = size(EEList,2);
    EEListOld = EEList;
    EEList = zeros(numEn,lenStre+datastep_add);
    EEList(:,1:lenStre) = EEListOld;
    clear('EEListOld')
    IeffListOld = IeffList;
    IeffList = zeros(5,numel,lenStre+datastep_add);
    IeffList(:,:,1:lenStre) = IeffListOld;
    clear('IeffListOld')
end

% Implicit error
% if impliciterr
%     if ~exist('IEHList','var')
%         IEHList = stepmax;
%     end
% %         IEList = zeros(***,length(EEHList)); % implicit error is
% %         specially implemented
%     % start with first step that the user desires quantities
%     stepIE = 1;
% end

% Plasticity
if PHist == 1
    PlasListOld = PlasList;
    PlasList = zeros(12,numel*nen,datastep_new);
    PlasList(:,1:numel*nen,1:datastep) = PlasListOld;
    clear('PlasListOld')
end

% Energy
if EHist == 1
    EnerListOld = EnerList;
    EnerList = zeros(1,datastep_new);
    EnerList(:,1:datastep) = EnerListOld;
    clear('EnerListOld')
end

% macro stresses
if SSHist == 1
    SSListOld = SSList;
    SSList = zeros(13,datastep_new);
    SSList(:,1:datastep) = SSListOld;
    clear('SSListOld')
end

% Interface quantities
if IFHist == 1
    InterQuantOld = InterQuant;
    InterQuant = zeros(numIQ,numnpI,datastep_new);
    InterQuant(:,1:datastep) = InterQuantOld;
    clear('InterQuantOld')
end

datastep = datastep_new;
