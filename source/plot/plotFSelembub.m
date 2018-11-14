function plotFSelembub(Coordinates, Contour, ndm, nen, ModelID, subp, numsubp, plotname,nelPn,zcont3D)
% plotModel(Coordinates, PatchTable, PatchIndex, PatchConn, U, V, W,
% numel, ModelID, geomref)
%
% Tim Truster
% CEE Graduate Student
% UIUC
% 07/16/2013
%
% Plot the edge bubble functions for fine scales
%close all

% Default arguments
if nargin == 5
    subp = 1;
    numsubp = 1;
    plotname = '';
    nelPn = [3 4 6 9 0];
    zcont3D = 0;
elseif nargin == 6
    disp('must specify both subplot # and number of subplots')
    return
elseif nargin == 7
    plotname = '';
    nelPn = [3 4 6 9 0];
    zcont3D = 0;
elseif nargin == 8
    nelPn = [3 4 6 9 0];
    zcont3D = 0;
elseif nargin == 9
    zcont3D = 0;
end

%Select and clear figure
figure(ModelID)
% set(ModelID,'Position',[360   425   560   273])
% set(ModelID,'Position',[360   278   560   420]) %old HP
set(ModelID,'Position',[403   246   560   420])
subplot(numsubp,1,subp)

colormap jet
% facecolors = [1 1 0     %yellow  U1 1
%               0 0.8 0   %green   U2 2
%               1 0 1     %magenta V1 3
%               0.8 0 0   %red     V2 4
%               0 1 1     %cyan    W1 5
%               0 0 1];   %blue    W2 6

%Plot patch faces from ModelID

nelP3 = nelPn(1);
nelP4 = nelPn(2);
nelP6 = nelPn(3);
nelP9 = nelPn(4);
contP = nelPn(5);

if ndm == 3
   
    
else
    
%     for elem = 1:numel
    %     pause
    
        %Determine patch size parameters
        nel = nen;
        
        if nel == 3
            nelP = nelP3;
        elseif nel == 4
            nelP = nelP4;
        elseif nel == 6
            nelP = nelP6;
        elseif nel == 9
            nelP = nelP9;
        end
        if contP == 1
            nel = nelP;
        end

            %Plot W Face
%             PatchPoints = zeros(nel,2);
%             ContourValues = zeros(nel,1);
%             
%             %Extract face points
%             for j = 1:nel
%                 node = ix(elem,j);
%                 for l = 1:2
                    PatchPoints = Coordinates;
%                 end %l
                ContourValues = Contour;
%             end %j
            
            %Plot Face
            %Plot mesh-less face
%             if geomref == 'g'
        if nel == 3 || nel == 6
%             len = 3;
            len = 5;
%             knots = [0 0.50 1.0 0 2/2/4 1/2 0 0 0
%                      0 0 0 0.5  0.5  0.5 1.0  1.0  1.0]';
            knots = [0 0.25 0.50 0.75 1.0 0 3/4/4 6/4/4 9/4/4 12/4/4 0 1/2/4 2/2/4 3/2/4 1/2 0 1/4/4 2/4/4 3/4/4 1/4 0 0 0 0 0
                     0 0 0 0 0 0.25 0.25 0.25 0.25 0.25 0.5  0.5  0.5  0.5  0.5 0.75 0.75 0.75 0.75 0.75 1.0  1.0  1.0  1.0  1.0]';

            plot2DLagransurfebub(PatchPoints, ContourValues, nel, knots, len, 'none', 0,zcont3D)

        elseif nel == 4 || nel == 9
            len = 5;
            knots = [-1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0
                     -1.0 -1.0 -1.0 -1.0 -1.0 -0.5 -0.5 -0.5 -0.5 -0.5 -0.0 -0.0 -0.0 -0.0 -0.0  0.5  0.5  0.5  0.5  0.5  1.0  1.0  1.0  1.0  1.0]';

            plot2DLagransurfebub(PatchPoints, ContourValues, nel, knots, len, 'none', 0,zcont3D)
%             end

        end

%     end
%     colorbar('EastOutside')
%     title(plotname,'FontSize',14)
    shading interp
    xlabel x
    ylabel y
%     zlabel z

% details for figures for Masters Thesis
% axis square
% axis([-0.1 10.1 -2 2])
axis equal
% axis([0 10 -3 3])
axis off
% axis([-0.1  40.1  -8.8  22.8])
colorbar('EastOutside','FontSize',13,'FontName','Arial','LineWidth',1)
% colorbar('EastOutside','FontSize',13,'FontName','Arial','LineWidth',1,'YTick',1e4*[-5 -2.5 0 2.5 5],'YTickLabel',[-5 -2.5 0 2.5 5])
    
end