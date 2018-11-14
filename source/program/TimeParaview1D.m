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


% This version assumes consecutive node numbering; the line elements are
% converted into a quadrilateral for plotting to make it easier to see

%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN CONTROL BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steps to output: starting step, increment to step, and last step
stepstart = 0; % can be zero if the initial '0' arrays are defined
stepinc = 1;
stepstop = 100;

% Nodes to output
nodestart = 1;
nodeend = numnp;

% Elements to output
elemstart = 1;
elemend = numel;

% Number of nodes per element; assumed constant for all elements elemstart 
% to elemend
nenPV = nen; % nen to use for Paraview file; change from nen when DG elements are used
thickness = 1; % thickness of bar for plotting

strflag = 0; % flag for stress postprocessing
serflag = 1; % flag for elemental stress postprocessing
errflag = 0; % flag for local-explicit error processing
preflag = 0; % flag for pressure
disflag = 1; % flag for displacement
velflag = 1; % flag for velocity
accflag = 1; % flag for acceleration
PMflag = 1; % flag for printing particle/matrix colors

dfoldername = 'ParaOut'; %default name of folder for files to be placed
userdefaultdir = 'C:\Users\Truster\Documents\Research\IO_Files_and_Codes\IO_Files-Tim\output'; %list your default output root directory here
foldername = uigetdir(userdefaultdir,'Select folder to write VTK files');
filename = 'out'; %prefix filename (e.g. 'out' -> out01.vtk)
header = 'Elasticity'; %text at top of VTK files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END CONTROL BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ndm == 1
nelP3 = 2;
nelP4 = 3;
nelP6 = 4;
nelP9 = 5;
nelS3 = 2;
nelS4 = 3;
nelS6 = 4;
nelS9 = 5;
else
    error('this version is only for 1D analysis')
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

numnpe = nodeend - nodestart + 1;
numele = elemend - elemstart + 1;
NodeList1 = NodesOnElement(elemstart:elemend,[2 1]);
NodeList2 = NodeList1 + numnp;
NodesOnElement1 = [NodesOnElement(elemstart:elemend,1:2) NodeList2];
RegionOnElement1 = RegionOnElement(elemstart:elemend);
Coordinates1 = [Coordinates -thickness/2*ones(numnp,1)
             Coordinates  thickness/2*ones(numnp,1)];
nodeend = 2*numnpe;
numnpe = nodeend - nodestart + 1;

%% Print results to VTK file: loop over steps
for stepi = stepstart:stepinc:stepstop
    
    if stepstop < 10000
        if stepi < 10
            parfilename = strcat(foldername,'/',filename,'000',num2str(stepi));
        elseif stepi < 100
            parfilename = strcat(foldername,'/',filename,'00',num2str(stepi));
        elseif stepi < 1000
            parfilename = strcat(foldername,'/',filename,'0',num2str(stepi));
        else
            parfilename = strcat(foldername,'/',filename,num2str(stepi));
        end
    end
    
    if stepi == 0
    NodeDispPV = DispList0(:,[1:numnp 1:numnp]);
    else
    NodeDispPV = squeeze(DispList(:,[1:numnp 1:numnp],stepi));
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
            StreList3 = zeros(npstr,2*numnpe);
        else
            if stepi > 0
                StreList3 = squeeze(StreList(1:npstr,[1:numnp 1:numnp],stepi));
            else
                StreList3 = zeros(npstr,2*numnpe);
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
            celltype = 9;
            nen2 = 4;
            elemtext = ' %i %i %i %i %i\n';
        otherwise
            error('only linear elements tested')

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
        fprintf(fid,'%16.7E %16.7E %16.7E\n',[Coordinates1(nodestart:nodeend,1:2)'; zeros(1,numnpe)]);
%     end
    end
%% Elements
    fprintf(fid,'\n');
    fprintf(fid,'CELLS %1$i %2$i\n',numele,(nen2+1)*numele);

%     for elem = 1:numel
        fprintf(fid,elemtext,[nen2*ones(numele,1) (NodesOnElement1(elemstart:elemend,1:nen2)-(nodestart)*ones(numele,nen2))]');
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
          end
%         end
end
%% Velocity
if velflag && exist('VeloList','var')
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Velocity float\n');

    if stepi == 0
    NodeDispPV = VeloList0(:,[1:numnp 1:numnp]);
    else
    NodeDispPV = squeeze(VeloList(:,[1:numnp 1:numnp],stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 1
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1,nodestart:nodeend); zeros(2,numnpe)]);
        end
end
%% Acceleration
if accflag && exist('AcceList','var')
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Acceleration float\n');

    if stepi == 0
    NodeDispPV = AcceList0(:,[1:numnp 1:numnp]);
    else
    NodeDispPV = squeeze(AcceList(:,[1:numnp 1:numnp],stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 1
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[NodeDispPV(1,nodestart:nodeend); zeros(2,numnpe)]);
        end
end
%% Stresses (xx,yy,xy)
%     if strflag == 1
%         
%         if ndm == 1
% 
%             
%         fprintf(fid,'\n');
%         fprintf(fid,'POINT_DATA   %i\n',numnpe);
%         fprintf(fid,'SCALARS Stress float\n');
%         fprintf(fid,'LOOKUP_TABLE default\n');
% 
%     if stepi == 0
%     NodeStrePV = zeros(1,2*numnp);
%     else
%     NodeStrePV = squeeze(StreList3(:,[1:numnp 1:numnp],stepi));
%     end
%     
%         fprintf(fid,'%16.7E\n',NodeStrePV);
%             
%         
%         end
%     end
    

        fprintf(fid,'\n');
        fprintf(fid,'CELL_DATA   %i\n',numele);
        
%% Elemental Stresses
    if serflag == 1 && exist('StreListE','var')
        
        if ndm == 1

            
        fprintf(fid,'\n');
%         fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS EStress float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

    if stepi == 0
    ElemStrePV = zeros(1,numele);
    else
    ElemStrePV = squeeze(StreListE(2,1:numele,stepi));
    end

        fprintf(fid,'%16.7E\n',ElemStrePV);
            
        
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
    
end
