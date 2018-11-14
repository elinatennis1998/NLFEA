% Tim Truster
% 04/14/2014
%
% Generate a FEAP input file from MATLAB mesh/problem data
% Assumes that the problem information (e.g. input file) have already been
% loaded to the Workspace.


filename = input('Specify filename for input file, excluding I: ','s');
filename = [filename '.txt'];
% Change file system designation based on computer type
if ispc
    directory = [pwd '\'];
elseif isunix
    directory = [pwd '/'];
else
    error('not sure what computer this is')
end


% Run minimum amount of FE data structure generation routines, to compute
% nodal loads
lamda = 1;

% Read major problem size parameters and solution options
pstart

% Initialize default parameters
pdefault

% Set pointers for allocation of mesh arrays
nie = ndf + 8;

% Allocate and zero arrays
ieFEAP = zeros(nie,nummat);
iedof = zeros(ndf,nen,nummat);
% ixFEAP = zeros(7,numel);
% ixFEAP(7,:) = RegionOnElement;

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


%% ADD SOMETHING HERE TO FIND OUT HOW MANY LOADS ARE ACTUALLY NONZERO
% Also treat np and regular loads/BC differently
numF = numnp;



% Set up element connectivity printing based on nen
if nen <= 13
elemtext = repmat(' %i',[1,nen]);
elemtext = ['%i 0 %i',elemtext,' \n'];
elseif nen <= 13+16
elemtext1 = repmat(' %i',[1,13]);
elemtext2 = repmat(' %i',[1,nen-13]);
elemtext = ['%i 0 %i',elemtext1,' \n',elemtext2,' \n'];
elseif nen <= 13+16+16
elemtext1 = repmat(' %i',[1,13]);
elemtext2 = repmat(' %i',[1,16]);
elemtext3 = repmat(' %i',[1,nen-13-16]);
elemtext = ['%i 0 %i',elemtext1,' \n',elemtext2,' \n',elemtext3,' \n'];
else
    error('too many nodes per element');
end

% Boundary and Force formats
boun2 = '%i,0,%i,%i\n';
boun3 = '%i,0,%i,%i,%i\n';
boun4 = '%i,0,%i,%i,%i,%i\n';
disp2 = '%i,0,%1.13f,%1.13f\n';
disp3 = '%i,0,%1.13f,%1.13f,%1.13f\n';
disp4 = '%i,0,%1.13f,%1.13f,%1.13f,%1.13f\n';
forc2 = '%i,0,%1.13f,%1.13f\n';
forc3 = '%i,0,%1.13f,%1.13f,%1.13f\n';
forc4 = '%i,0,%1.13f,%1.13f,%1.13f,%1.13f\n';
if ndf == 2
    boun = boun2;
    displ = disp2;
    forc = forc2;
elseif ndf == 3
    boun = boun3;
    displ = disp3;
    forc = forc3;
elseif ndf == 4
    boun = boun4;
    displ = disp4;
    forc = forc4;
end

NodeBC = sortrows(NodeBC);

fid = fopen([directory 'I' filename],'wt');
if numBC > 20
fidB = fopen([directory 'B' filename],'wt');
else
fidB = fid;
end
if numnp > 100
fidC = fopen([directory 'C' filename],'wt');
else
fidC = fid;
end
if numel > 100
fidE = fopen([directory 'E' filename],'wt');
else
fidE = fid;
end
if numComp > 100
fidT = fopen([directory 'T' filename],'wt');
else
fidT = fid;
end
if numF > 100
fidV = fopen([directory 'V' filename],'wt');
else
fidV = fid;
end
fidM = fopen([directory 'Mat' filename],'wt');
fidD = fopen([directory 'D' filename],'wt');
fidU = fopen([directory 'U' filename],'wt');
fprintf(fid,'feap * Name Here *\n');
fprintf(fid,'%i,%i,%i,%i,%i,%i\n',ProbType);
fprintf(fid,'\n');
fprintf(fid,'NOPRint\n');
fprintf(fid,'\n');
fprintf(fid,'COORds\n');
if fidC ~= fid
fprintf(fid,'INCLude %s\n',['C' filename]);
end
if ndm == 2
fprintf(fidC,'%i %i %16.7E %16.7E\n',[(1:numnp); zeros(1,numnp); Coordinates']);
elseif ndm == 3
fprintf(fidC,'%i %i %16.7E %16.7E %16.7E\n',[(1:numnp); zeros(1,numnp); Coordinates']);
else
fprintf(fidC,'%i %i %16.7E\n',[(1:numnp); zeros(2,numnp); Coordinates']);
end
fprintf(fid,'\n');
fprintf(fid,'ELEMents\n');
if fidE ~= fid
fprintf(fid,'INCLude %s\n',['E' filename]);
end
fprintf(fidE,elemtext,[(1:numel); RegionOnElement(1:numel)'; NodesOnElement(1:numel,1:nen)']);

% REMOVE THIS SOON
if numSI > 0
fidN = fopen([directory 'N' filename],'wt');
fprintf(fid,'\n');
fprintf(fid,'INCLude %s\n',['N' filename]);
fprintf(fidN,'SURI,%i\n',numSI);
fprintf(fidN,' %i %i %i %i %i %i %i %i %i %i\n',[SurfacesI(:,1:8) ones(numSI,1)*MateT(nummat,3:4)]');
end

%% Boundary Conditions

NodeBC2 = zeros(numnp,ndf);
NodeBC3 = zeros(numnp,ndf);

dflag = 0;

for bc = 1:numBC
    
    node = NodeBC(bc,1);
    dof = NodeBC(bc,2);
    NodeBC2(node,dof) = 1;
    if NodeBC(bc,3) ~= 0
        dflag = 1;
        NodeBC3(node,dof) = NodeBC(bc,3);
    end
    
end

FNodeList = zeros(numnp,ndf);
% Store groups of nodes for either displacement or force: proportional
% loading

% Print header
fprintf(fid,'\n');
fprintf(fid,'BOUN\n');
if fidB ~= fid
fprintf(fid,'INCLude %s\n',['B' filename]);
end
for node = 1:numnp
    
    if sum(NodeBC2(node,:)) > 0
        fprintf(fidB,boun,node,NodeBC2(node,1:ndf));
    end
    
end

fprintf(fid,'\n');

if dflag == 1
    
    % Print header
    fprintf(fid,'DISP\n');
    fprintf(fid,'INCLude %s\n',['D' filename]);
    
    for node = 1:numnp
    
        if sum(abs(NodeBC3(node,:))) > 0
            fprintf(fidD,displ,node,NodeBC3(node,1:ndf));
            for i = 1:ndf
                if NodeBC3(node,i) ~= 0
                    FNodeList(node,i) = 1;
                end
            end
        end

    end

    fprintf(fid,'\n');
    
end

fprintf(fid,'\n');

% Proportional load groups: DISP = 1, numSL = 2, numSLnp = 3
fprintf(fid,'! For displacements in NodeBC\n');
fprintf(fid,'FPROp\n');
if fidU ~= fid
fprintf(fid,'INCLude %s\n',['U' filename]);
end

for node = 1:numnp
    
    if sum(abs(FNodeList(node,:))) > 0
        if ndf == 3
        fprintf(fidU,'%i,0,%i,%i,%i\n',node,FNodeList(node,:));
        else
            disp('bad ndf')
            return
        end
    end

end
FNodeList = zeros(numnp,ndf);

% Print nodal equivalent forces (proportional)
% Header
fprintf(fid,'! All forces except numSLnp\n');
fprintf(fid,'FORC\n');

NodeLoad2 = zeros(numnp,ndf);

for node = 1:numnp

    for dir = 1:ndf
        gDOF = NDOFT(node,dir);
        if gDOF <= neq && gDOF > 0
            NodeLoad2(node,dir) = FcU(gDOF); %%%%%%% Fc1 or FcU
        end
    end

end

for node = 1:numnp
    
    if sum(abs(NodeLoad2(node,:))) > 0
        fprintf(fid,forc,node,NodeLoad2(node,1:ndf));
        for i = 1:ndf
            if NodeLoad2(node,i) ~= 0
                FNodeList(node,i) = 2;
            end
        end
    end

end

% fprintf(fid,'\n');
% 
% % Print nodal equivalent forces (nonproportional)
% % Header
% fprintf(fid,'! Forces from numSLnp\n');
% fprintf(fid,'FORC\n');
% 
% NodeLoad2 = zeros(numnp,ndf);
% 
% for node = 1:numnp
% 
%     for dir = 1:ndf
%         gDOF = NDOFT(node,dir);
%         if gDOF <= neq && gDOF > 0
%             NodeLoad2(node,dir) = Fc1np(gDOF);
%         for i = 1:ndf
%             if NodeLoad2(node,i) ~= 0
%                 FNodeList(node,i) = 3;
%             end
%         end
%         end
%     end
% 
% end
% 
% for node = 1:numnp
%     
%     if sum(abs(NodeLoad2(node,:))) > 0
%         fprintf(fid,forc,node,NodeLoad2(node,1:ndf));
%     end
% 
% end

% fprintf(fid,'\n');
% 
% % Proportional load groups: DISP = 1, numSL = 2, numSLnp = 3
% fprintf(fid,'FPROp\n');
% if fidV ~= fid
% fprintf(fid,'INCLude %s\n',['V' filename]);
% end
% 
% for node = 1:numnp
%     
%     if sum(abs(FNodeList(node,:))) > 0
%         if ndf == 3
%         fprintf(fidV,'%i,0,%i,%i,%i\n',node,FNodeList(node,:));
%         else
%             disp('bad ndf')
%             return
%         end
%     end
% 
% end

fprintf(fid,'\n');
fprintf(fid,'INCLude %s\n',['Mat' filename]);

%% Material and Solution Data
% Likely needs changed for each different problem
doflist = '';
for j = 1:ndf
    doflist = [doflist ' %i'];
end
proplist = '';
for j = 1:size(MateT,2)
    proplist = [proplist ' %5.6e'];
end
for ma = 1:nummat
    fprintf(fidM,'MATErial %i\n',ma);
    if size(MatTypeTable,1) > 3
        text = ['  USER %i %i ' doflist ' nonlinear %i\n'];
        fprintf(fidM,text,MatTypeTable([2 1 4:ndf+3 3],ma));
    else
        fprintf(fidM,'  USER %i %i nonlinear %i\n',MatTypeTable([2 1 3],ma));
    end
    text = ['  ' proplist '\n'];
    fprintf(fidM,text,MateT(ma,:));
    fprintf(fidM,'\n');
    
end

fprintf(fid,'\n');
fprintf(fid,'PRINt\n');
fprintf(fid,'END\n');

%% Linking nodes

if numComp > 0
    
fprintf(fid,'LINK\n');
if fidT ~= fid
fprintf(fid,'INCLude %s\n',['T' filename]);
end

if ndf == 2
fprintf(fidT,'%i,%i,,,%i,%i\n',NodeComp');
else
fprintf(fidT,'%i,%i,,,%i,%i,%i\n',NodeComp(:,[2 1 3 4 5])');
end

end

fprintf(fid,'\n');
fprintf(fid,'INTE\n');
fprintf(fid,'\n');
fprintf(fid,'Batch\n');
fprintf(fid,' opti\n');
fprintf(fid,'end\n');

fprintf(fid,'STOP\n');

fclose('all');
