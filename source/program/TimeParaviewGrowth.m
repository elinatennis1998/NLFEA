%  11/04/2012
% Revised 08/21/2013
% Renamed 09/03/2017
% Tim Truster

% Function to output time history to vtk file for Paraview viewing; for
% mixture and growth modeling

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
stepstart = 1; % can be zero if the initial '0' arrays are defined
stepinc = 1;
stepstop = 1;

% Nodes to output
nodestart = 1;
nodeend = numnp;

% Elements to output
elemstart = 1;
elemend = numel;

strflag = 1; % flag for stress postprocessing
st2flag = 0; % flag for stress vectors/scalars
errflag = 0; % flag for local-explicit error processing
preflag = 1; % flag for pressure
disflag = 1; % flag for displacement
velflag = 0; % flag for velocity
denflag = 1; % flag for density
accflag = 0; % flag for acceleration
PMflag = 1; % flag for printing particle/matrix colors

dfoldername = 'ParaOut'; %name of folder for files to be placed
userdefaultdir = 'E:\Users\trustetj\Documents\Research\IO_Files-Tim\output'; %list your default output root directory here
foldername = uigetdir(userdefaultdir,'Select folder to write VTK files');
filename = 'out'; %prefix filename (e.g. 'out' -> out01.vtk)
header = 'Elasticity'; %text at top of VTK files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END CONTROL BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numnpe = nodeend - nodestart + 1;
numele = elemend - elemstart + 1;
if ndm == 2
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

if isnumeric(foldername)
    foldername = dfoldername;
end
if exist(foldername,'dir') == 0
    mkdir(foldername);
end


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
    Node_U_V = DispList0;
    else
    Node_U_V = squeeze(DispList(:,:,stepi));
    end
    
    if strflag == 1
        if(exist('StreList','var')==0)
            if ndm == 2
                numstr = 7; % this value can be changed
            elseif ndm == 3
                numstr = 11; % this value can be changed
            end
        StreList2 = zeros(numnp,numstr);
        Eareas = zeros(numnp,1);

        isw = 25;
        FormS2

        StreList3 = (StreList2./(Eareas*ones(1,numstr)))';
        else
            if stepi > 0
        StreList3 = squeeze(StreList(1:numstr,:,stepi));
            else
        StreList3 = zeros(numstr,numnpe);
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
    switch nen

        case 3
            celltype = 5;
            nen2 = nen;
            elemtext = ' %i %i %i %i\n';
        case 4
            if ndm == 2
                celltype = 9;
            else
                celltype = 10;
            end
            nen2 = nen;
            elemtext = ' %i %i %i %i %i\n';
        case 6
            celltype = 22;
            nen2 = nen;
            elemtext = ' %i %i %i %i %i %i %i\n';
        case 8
            celltype = 12;
            nen2 = nen;
            elemtext = ' %i %i %i %i %i %i %i %i %i\n';
        case 9
            celltype = 23;
            nen2 = 8;
            elemtext = ' %i %i %i %i %i %i %i %i %i\n';
        case 10
            celltype = 24;
            nen2 = nen;
            elemtext = ' %i %i %i %i %i %i %i %i %i %i %i\n';
        case 27
            celltype = 25;
            nen2 = 20;
            elemtext = ' %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n';

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
%     for node = 1:numnp
        fprintf(fid,'%16.7E %16.7E %16.7E\n',[Coordinates(nodestart:nodeend,:)'; zeros(1,numnpe)]);
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
        fprintf(fid,'%16.7E\n',Node_U_V(ndf,nodestart:nodeend)');
%     end
end
%% Displacement
if disflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Displacement float\n');

%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
          if ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[Node_U_V(1:2,nodestart:nodeend); zeros(1,numnpe)]);
          else
            fprintf(fid,'%16.7E %16.7E %16.7E\n',Node_U_V(1:3,nodestart:nodeend));
          end
%         end
end
%% Velocity
if velflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Velocity float\n');

    if stepi == 0
    Node_U_V = VeloList0(:,:,stepi);
    else
    Node_U_V = squeeze(VeloList(:,:,stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[Node_U_V(1:2,nodestart:nodeend,stepi); zeros(1,numnpe)]);
        elseif ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',Node_U_V(1:3,nodestart:nodeend,stepi));
        end
end
%% Acceleration
if accflag
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS Acceleration float\n');

    if stepi == 0
    Node_U_V = AcceList0(:,:,stepi);
    else
    Node_U_V = squeeze(AcceList(:,:,stepi));
    end
%         for node = 1:numnp
    %         fprintf(fid,'%1$16.7E %2$16.7E %3$16.7E\n',ScaVecTen(1,node),ScaVecTen(2,node),ScaVecTen(3,node));
        if ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',[Node_U_V(1:2,nodestart:nodeend,stepi); zeros(1,numnpe)]);
        elseif ndm == 2
            fprintf(fid,'%16.7E %16.7E %16.7E\n',Node_U_V(1:3,nodestart:nodeend,stepi));
        end
end
%% Stresses (xx,yy,xy)
    if strflag == 1
        if ndm == 3
        
        if st2flag
            
            if numstr > 6
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS VonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(7,nodestart:nodeend));
        
            if numstr > 7
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS StressVec float\n');

            fprintf(fid,'%16.7E %16.7E %16.7E\n',StreList3(8:10,nodestart:nodeend));
            
            if numstr > 10
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
            if numstr > 3
    fprintf(fid,'\n');
%     fprintf(fid,'POINT_DATA   %i\n',numnp);
    fprintf(fid,'SCALARS VonMises float\n');
    fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',StreList3(4,nodestart:nodeend));
            if numstr > 4
        fprintf(fid,'\n');
        fprintf(fid,'VECTORS StressVec float\n');

            fprintf(fid,'%16.7E %16.7E 0.0000000E+000\n',StreList3(5:6,nodestart:nodeend));
            if numstr > 6
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
            
        end
    end

%% H1 local explicit error
    if errflag == 1

        fprintf(fid,'\n');
        fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS LocalExpError float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',sqrt(Ieffvals(:,2))');
    end


%% Colors for particle/matrix or material flag identifier
    if PMflag == 1

        fprintf(fid,'\n');
        fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS ParMat float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');

        fprintf(fid,'%16.7E\n',RegionOnElement(elemstart:elemend)');
    end

%% Density
    if denflag == 1 % Only works for NL_Elem2_3dMG.m

        fprintf(fid,'\n');
%         fprintf(fid,'CELL_DATA   %i\n',numele);
        fprintf(fid,'SCALARS Density float\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        
        Densities = hrmat(elemstart:elemend*16,stepi);
        Densities = reshape(Densities,16,numele);
        Densities = sum(Densities(1:8,:),1)'/8;
        Densities = (1 + Densities) * rho0;

        fprintf(fid,'%16.7E\n',Densities);
    end

    fclose(fid);
    
end
