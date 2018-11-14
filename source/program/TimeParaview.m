%  11/04/2012
% Revised 08/21/2013
% Tim Truster

% Function to output time history to vtk file for Paraview viewing

% README:
% 0. Specify in the variable "userdefaultdir" the path to the root directory
%    that you plan to put all of your output folders
% 1. Execute file with F5 or play button
% 2. A window opens to specify the output folder for the current VTK files;
%    to create a new folder within the "userdefaultdir" directory quickly,
%    follow the steps below
% 3. Click the "Make New Folder" button in the window, which will create a
%    new directory in the current folder (currently "userdefaultdir")
% 4. Press "Enter" key to confirm folder name
% 5. Press "Enter" key to select the new folder
% 6. Script then executes and writes the VTK files


%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN CONTROL BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steps to output: starting step, increment to step, and last step
stepstart = 0; % can be zero if the initial '0' arrays are defined
stepinc = 1;
stepstop =3000;100;500;1000;

startp = stepstart; % first step # for VTK file
incp = 1; % first step # for VTK file
stopp = stepstop/1; % first step # for VTK file

% Nodes to output
nodestart = 1;
nodeend = numnp;

% Elements to output
elemstart = 1;
elemend = numel-numSI;

% Number of nodes per element; assumed constant for all elements elemstart 
% to elemend
nenPV = nen/2;%/2;4;10; % nen to use for Paraview file; change from nen when DG elements are used

% Couplers to output
coupstart = numel-numSI+1;
coupend = numel;

strflag = 1; % flag for stress postprocessing
st2flag = 1; % flag for stress vectors/scalars
serflag = 0; % flag for elemental stress postprocessing
se2flag = 0; % flag for elemental stress vectors/scalars
errflag = 0; % flag for local-explicit error processing
preflag = 0; % flag for pressure
disflag = 1; % flag for displacement
velflag = 0; % flag for velocity
accflag = 0; % flag for acceleration
PMflag = 1; % flag for printing particle/matrix colors

CZDG_flag = 1;2; % flag for CZ=1 or DG=2
cflag = 1;0; % flag for couplers
cseflag = 1;0; % flag for coupler element stresses

dfoldername = 'ParaOut'; %default name of folder for files to be placed
userdefaultdir = 'C:\Users\Truster\Documents\Research\IO_Files_and_Codes\IO_Files-Tim\output'; %list your default output root directory here
foldername = uigetdir(userdefaultdir,'Select folder to write VTK files');
filename = 'out'; %prefix filename (e.g. 'out' -> out01.vtk)
header = 'Elasticity'; %text at top of VTK files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END CONTROL BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numnpe = nodeend - nodestart + 1;
numele = elemend - elemstart + 1;
numcoup = coupend - coupstart + 1;
if ndm == 1
nelP3 = 2;
nelP4 = 3;
nelP6 = 4;
nelP9 = 5;
nelS3 = 2;
nelS4 = 3;
nelS6 = 4;
nelS9 = 5;
elseif ndm == 2
nelP3 = 3;
nelP4 = 4;
nelP6 = 6;
nelP9 = 4;
nelS3 = 3;
nelS4 = 4;
nelS6 = 3;
nelS9 = 4;
else
nelP3 = 4;
nelP4 = 8;
nelP6 = 10;
nelP9 = 27;
nelS3 = 4;
nelS4 = 8;
nelS6 = 4;
nelS9 = 8;
end

% check if foldername is a number, meaning the user canceled and didn't
% pick a folder
if isnumeric(foldername)
    foldername = dfoldername;
end
% create the folder if it does not exist
if exist(foldername,'dir') == 0
    mkdir(foldername);
end

% ensure that the output step requested does not exceed the existing
% analysis output
if exist('stepmax','var')
    stepstop = min(stepstop,stepmax);
end
% make sure the initial step arrays do exist
if ~exist('DispList0','var')
    stepstart = max(stepstart,1);
end

%% Print results to VTK file: loop over steps
paravsteps = startp:incp:stopp;
stept = 0;
for stepi = stepstart:stepinc:stepstop
    
    stept = stept+1;
    stepp = paravsteps(stept);
    
    if stepstop < 10000
        if stepi < 10
            parfilename = strcat(foldername,'/',filename,'000',num2str(stepp));
        elseif stepi < 100
            parfilename = strcat(foldername,'/',filename,'00',num2str(stepp));
        elseif stepi < 1000
            parfilename = strcat(foldername,'/',filename,'0',num2str(stepp));
        else
            parfilename = strcat(foldername,'/',filename,num2str(stepp));
        end
    end
    
    if stepi == 0
    NodeDispPV = DispList0;
    else
    NodeDispPV = squeeze(DispList(:,:,stepi));
    end
    
    if strflag == 1
        if(exist('StreList','var')==0)
%             if ndm == 2
%                 npstr = 7; % this value can be changed
%             elseif ndm == 3
%                 npstr = 11; % this value can be changed
%             elseif ndm == 1
%                 npstr = 1; % this value can be changed
%             end
%         StreList2 = zeros(numnp,npstr);
%         Eareas = zeros(numnp,1);
% 
%         isw = 25;
%         FormS2
% 
%         StreList3 = (StreList2./(Eareas*ones(1,npstr)))';
            disp('warning, StreList is not defined from analysis; zeros are used')
            StreList3 = zeros(npstr,numnpe);
        else
            if stepi > 0
                StreList3 = squeeze(StreList(1:npstr,:,stepi));
            else
                StreList3 = zeros(npstr,numnpe);
            end
        end
    end
    
    if serflag == 1
        if(exist('StreListE','var')==0)
%             if ndm == 2
%                 npstr = 7; % this value can be changed
%             elseif ndm == 3
%                 npstr = 11; % this value can be changed
%             elseif ndm == 1
%                 npstr = 1; % this value can be changed
%             end
%         StreList2 = zeros(numnp,npstr);
%         Eareas = zeros(numnp,1);
% 
%         isw = 25;
%         FormS2
% 
%         StreList3 = (StreList2./(Eareas*ones(1,npstr)))';
            disp('warning, StreListE is not defined from analysis; zeros are used')
            StreList4 = zeros(nestr,numele);
        else
            if stepi > 0
                StreList4 = squeeze(StreListE(1:nestr,:,stepi));
            else
                StreList4 = zeros(nestr,numele);
            end
        end
    end
%     if errflag == 1
%         bubblevals = zeros(numel,3);
%         Ieffvals = zeros(numel,3);
% 
%         isw = 11;
%         FormE
% 
%     end
    

    % Set cell-type and number of nodes on elements
    switch nenPV

        case 2
            celltype = 3;
            nen2 = nenPV;
            elemtext = ' %i %i %i\n';
        case 3
            celltype = 5;
            nen2 = nenPV;
            elemtext = ' %i %i %i %i\n';
            intertype = 3;
            nenI = 2;
            couptext = ' %i %i %i\n';
            couporder = [1 2];
        case 4
            if ndm == 2
                celltype = 9;
                intertype = 3;
                nenI = 2;
                couptext = ' %i %i %i\n';
                couporder = [1 2];
            elseif ndm == 3
                celltype = 10;
                intertype = 5;
                nenI = 3;
                couptext = ' %i %i %i %i\n';
                couporder = [1 2 3];
            end
            nen2 = nenPV;
            elemtext = ' %i %i %i %i %i\n';
        case 6
            celltype = 22;
            nen2 = nenPV;
            elemtext = ' %i %i %i %i %i %i %i\n';
            intertype = 21;
            nenI = 3;
            couptext = ' %i %i %i %i\n';
            couporder = [1 2 4];
        case 8
            if ndm == 3
                celltype = 12;
                intertype = 9;
                nenI = 4;
                couptext = ' %i %i %i %i %i\n';
                couporder = [1 2 3 4];
            elseif ndm == 2
                celltype = 23;
                intertype = 21;
                nenI = 3;
                couptext = ' %i %i %i %i\n';
                couporder = [1 2 5];
            end
            nen2 = nenPV;
            elemtext = ' %i %i %i %i %i %i %i %i %i\n';
        case 9
            celltype = 23;
            nen2 = 8;
            elemtext = ' %i %i %i %i %i %i %i %i %i\n';
            intertype = 21;
            nenI = 3;
            couptext = ' %i %i %i %i\n';
            couporder = [1 2 5];
        case 10
            celltype = 24;
            nen2 = nenPV;
            elemtext = ' %i %i %i %i %i %i %i %i %i %i %i\n';
            intertype = 22;
            nenI = 6;
            couptext = ' %i %i %i %i %i %i %i\n';
            if CZDG_flag == 1
            couporder = [1 2 3 4 5 6];
            else
            couporder = [1 2 3 5 6 7];
            end
        case 20
            celltype = 25;
            nen2 = nenPV;
            elemtext = ' %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n';
            intertype = 23;
            nenI = 8;
            couptext = ' %i %i %i %i %i %i %i %i %i\n';
            couporder = [1 2 3 4 9 10 11 12];
        case 27
            celltype = 25;
            nen2 = 20;
            elemtext = ' %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n';
            intertype = 23;
            nenI = 8;
            couptext = ' %i %i %i %i %i %i %i %i %i\n';
            couporder = [1 2 3 4 9 10 11 12];
            
    end

    fid = fopen(strcat(parfilename,'.vtk'),'wt');

    fprintf(fid,'# vtk DataFile Version 1.0\n');
    fprintf(fid,'%s \n',header);
    fprintf(fid,'ASCII\n');
    fprintf(fid,'\n');

    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS   %i float\n',numnpe);
%% Coordinates
    if ndm == 1
%     for node = 1:numnp
        fprintf(fid,'%16.7E %16.7E %16.7E\n',[Coordinates(nodestart:nodeend,1)'; zeros(2,numnpe)]);
%     end
    elseif ndm == 2
%     for node = 1:numnp
        fprintf(fid,'%16.7E %16.7E %16.7E\n',[Coordinates(nodestart:nodeend,1:2)'; zeros(1,numnpe)]);
%     end
    else
%     for node = 1:numnp
        fprintf(fid,'%16.7E %16.7E %16.7E\n',Coordinates(nodestart:nodeend,:)');
%     end
    end
%% Elements
    fprintf(fid,'\n');
    fprintf(fid,'CELLS %1$i %2$i\n',numele,(nen2+1)*numele);

%     for elem = 1:numel
        fprintf(fid,elemtext,[nen2*ones(numele,1) (NodesOnElement(elemstart:elemend,1:nen2)-(nodestart)*ones(numele,nen2))]');
%         for node = 1:nen2-1
%             fprintf(fid,' %i',NodesOnElement(elem,node)-1);
%         end
%         fprintf(fid,' %i\n',NodesOnElement(elem,nen2)-1);
%     end

    fprintf(fid,'\n');
    fprintf(fid,'CELL_TYPES %1$i\n',numele);

%     for elem = 1:numel
        fprintf(fid,'%i\n',celltype*ones(numele,1));
%     end

        fprintf(fid,'\n');
        fprintf(fid,'POINT_DATA   %i\n',numnpe);

%% Pressure
if preflag
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

%     for node = 1:numnp
        fprintf(fid,'%16.7E\n',NodeDispPV(ndf,nodestart:nodeend)');
%     end
end
%% Displacement
if disflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Displacement float\n');

%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
          if ndm == 1
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1,nodestart:nodeend); zeros(2,numnpe)]);
          elseif ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1:2,nodestart:nodeend); zeros(1,numnpe)]);
          else
            fprintf(fid,'%16.7E %16.7E %16.7E\n',NodeDispPV(1:3,nodestart:nodeend));
          end
%         end
end
%% Velocity
if velflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Velocity float\n');

    if stepi == 0
    NodeDispPV = VeloList0;
    else
    NodeDispPV = squeeze(VeloList(:,:,stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 1
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1,nodestart:nodeend); zeros(2,numnpe)]);
        elseif ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1:2,nodestart:nodeend); zeros(1,numnpe)]);
        elseif ndm == 3
            fprintf(fid,'%16.7E %16.7E %16.7E\n',NodeDispPV(1:3,nodestart:nodeend));
        end
end
%% Acceleration
if accflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Acceleration float\n');

    if stepi == 0
    NodeDispPV = AcceList0;
    else
    NodeDispPV = squeeze(AcceList(:,:,stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 1
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1,nodestart:nodeend); zeros(2,numnpe)]);
        elseif ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1:2,nodestart:nodeend); zeros(1,numnpe)]);
        elseif ndm == 3
            fprintf(fid,'%16.7E %16.7E %16.7E\n',NodeDispPV(1:3,nodestart:nodeend));
        end
end
%% Stresses (xx,yy,xy)
    if strflag == 1
        
        if ndm == 3
        
            
        if st2flag
            
            if npstr > 6
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS VonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(7,nodestart:nodeend));
        
            if npstr > 7
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS StressVec float\n');

            fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList3(8:10,nodestart:nodeend));
            
            if npstr > 10
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS Hydro float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(11,nodestart:nodeend));
            end; end; end
        
        end

        fprintf(fid,'\n');
        fprintf(fid,'TENSORS Stress float\n');

        fprintf(fid,'%16.7E %16.7E %16.7E\n %16.7E %16.7E %16.7E\n %16.7E %16.7E %16.7E\n \n',StreList3([1 4 6 4 2 5 6 5 3],:));
        
        
        elseif ndm == 2
        
            
        if st2flag
            if npstr > 3
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS VonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(4,nodestart:nodeend));
            if npstr > 4
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS StressVec float\n');

            fprintf(fid,'%16.7E %16.7E 0.0000000E+000\n',StreList3(5:6,nodestart:nodeend));
            if npstr > 6
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS Hydro float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(7,nodestart:nodeend));
            
            end; end; end
        end

        fprintf(fid,'\n');
        fprintf(fid,'TENSORS Stress float\n');

        fprintf(fid,'%16.7E %16.7E 0.0000000E+000\n %16.7E %16.7E 0.0000000E+000\n 0.0000000E+000 0.0000000E+000 0.0000000E+000\n \n',StreList3([1 3 3 2],:));
        
        
        elseif ndm == 1

            
        fprintf(fid,'\n');
        fprintf(fid,'SCALARS Stress float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(1,:));
            
        
        end
    end

    
    
        fprintf(fid,'\n');
        fprintf(fid,'CELL_DATA   %i\n',numele);
        
%% Elemental Stresses
    if serflag == 1
        
        
        if ndm == 3
        
            
        if st2flag
            
            if nestr > 6
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS EVonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList4(7,elemstart:elemend));
        
            if nestr > 7
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS StressVec float\n');

            fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(8:10,elemstart:elemend));
            
            if nestr > 10
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS EHydro float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList4(11,elemstart:elemend));
            end; end; end
        
        end

        fprintf(fid,'\n');
        fprintf(fid,'TENSORS EStress float\n');

        fprintf(fid,'%16.7E %16.7E %16.7E\n %16.7E %16.7E %16.7E\n %16.7E %16.7E %16.7E\n \n',StreList4([1 4 6 4 2 5 6 5 3],:));
        
        
        elseif ndm == 2
        
            
        if st2flag
            if nestr > 3
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS EVonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList4(4,elemstart:elemend));
            if nestr > 4
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS EStressVec float\n');

            fprintf(fid,'%16.7E %16.7E 0.0000000E+000\n',StreList4(5:6,elemstart:elemend));
            if nestr > 6
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS EHydro float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList4(7,elemstart:elemend));
            
            end; end; end
        end

        fprintf(fid,'\n');
        fprintf(fid,'TENSORS EStress float\n');

        fprintf(fid,'%16.7E %16.7E 0.0000000E+000\n %16.7E %16.7E 0.0000000E+000\n 0.0000000E+000 0.0000000E+000 0.0000000E+000\n \n',StreList4([1 3 3 2],elemstart:elemend));
        
        
        elseif ndm == 1

            
        fprintf(fid,'\n');
        fprintf(fid,'SCALARS EStress float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList4(1,:));
            
        
        end
    end

%% H1 local explicit error
    if errflag == 1

        fprintf(fid,'\n');
%         fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS LocalExpError float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',sqrt(Ieffvals(:,2))');
    end


%% Colors for particle/matrix or material flag identifier
    if PMflag == 1

        fprintf(fid,'\n');
%         fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS ParMat float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',RegionOnElement(elemstart:elemend)');
    end

    fclose(fid);

%% Couplers
    if cflag
        
    if stepstop < 10000
        if stepi < 10
            parfilename = strcat(foldername,'/','c',filename,'000',num2str(stepp));
        elseif stepi < 100
            parfilename = strcat(foldername,'/','c',filename,'00',num2str(stepp));
        elseif stepi < 1000
            parfilename = strcat(foldername,'/','c',filename,'0',num2str(stepp));
        else
            parfilename = strcat(foldername,'/','c',filename,num2str(stepp));
        end
    end
    fid = fopen(strcat(parfilename,'.vtk'),'wt');

    fprintf(fid,'# vtk DataFile Version 1.0\n');
    fprintf(fid,'%s \n',header);
    fprintf(fid,'ASCII\n');
    fprintf(fid,'\n');

    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS   %i float\n',numnpe);
%% Coordinates
    if ndm == 2
        fprintf(fid,'%16.7E %16.7E %16.7E\n',[Coordinates(nodestart:nodeend,1:2)'; zeros(1,numnpe)]);
    else
        fprintf(fid,'%16.7E %16.7E %16.7E\n',Coordinates(nodestart:nodeend,:)');
    end
%% Couplers
    fprintf(fid,'\n');
    fprintf(fid,'CELLS %1$i %2$i\n',numcoup,(nenI+1)*numcoup);

    fprintf(fid,couptext,[nenI*ones(numcoup,1) (NodesOnElement(coupstart:coupend,couporder)-(nodestart)*ones(numcoup,nenI))]');

    fprintf(fid,'\n');
    fprintf(fid,'CELL_TYPES %1$i\n',numcoup);

    fprintf(fid,'%i\n',intertype*ones(numcoup,1));

    fprintf(fid,'\n');
    fprintf(fid,'CELL_DATA   %i\n',numcoup);

%% Colors for material flag identifier
    fprintf(fid,'\n');
    fprintf(fid,'SCALARS IRegion float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

    fprintf(fid,'%16.7E\n',RegionOnElement(coupstart:coupend)');
    
%% Coupler Stresses
    if cseflag == 1
        
        if(exist('StreListE','var')==0)
            disp('warning, StreListE is not defined from analysis; zeros are used')
            StreList4 = zeros(nestr,numcoup);
        elseif stepi == 0
            StreList4 = zeros(nestr,numcoup);
        else
            StreList4 = squeeze(StreListE(1:nestr,coupstart:coupend,stepi));
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'FIELD CoupEData 11\n');
        if ndm == 3
        fprintf(fid,'NormalVector 3 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(1:3,:));
        fprintf(fid,'StabilParam 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(4,:));
        fprintf(fid,'FEgap 3 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(5:7,:));
        fprintf(fid,'FEgap_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(8,:));
        fprintf(fid,'Zeta 3 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(9:11,:));
        fprintf(fid,'Zeta_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(12,:));
        fprintf(fid,'TotalFlux 3 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(13:15,:));
        fprintf(fid,'TF_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(16,:));
        fprintf(fid,'FEflux 3 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList4(17:19,:));
        fprintf(fid,'FEF_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(20,:));
        fprintf(fid,'Zeta_max 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(21,:));
        elseif ndm == 2
        fprintf(fid,'NormalVector 2 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E\n',StreList4(1:2,:));
        fprintf(fid,'StabilParam 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(3,:));
        fprintf(fid,'FEgap 2 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E\n',StreList4(4:5,:));
        fprintf(fid,'FEgap_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(6,:));
        fprintf(fid,'Zeta 2 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E\n',StreList4(7:8,:));
        fprintf(fid,'Zeta_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(9,:));
        fprintf(fid,'TotalFlux 2 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E\n',StreList4(10:11,:));
        fprintf(fid,'TF_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(12,:));
        fprintf(fid,'FEflux 2 %i float \n',numcoup);
        fprintf(fid,'%16.7E %16.7E\n',StreList4(13:14,:));
        fprintf(fid,'FEF_n 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(15,:));
        fprintf(fid,'Zeta_max 1 %i float \n',numcoup);
        fprintf(fid,'%16.7E\n',StreList4(16,:));
        end

    end
        
    fclose(fid);

    end
    
end
