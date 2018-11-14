function plotInterCont2(Coordinates, Contour, NodesOnElement, SurfacesI, numSI, PlotID, elemlist, MNE_ID, zcont3D, varargin)
% Tim Truster
% 04/01/2015
% 
% Plot of interface quantities for 2d meshes, Contour generated e.g. by
% elemental stresses from isw=26; see NL_Elem21_2d_2.
%
% Default usage (for copy/paste):
% plotInterCont2(Coordinates, squeeze(StreListE(13,:,5))', NodesOnElement, SurfacesI, numSI, 1, (1:numSI)', [1 0 0], 0)
%
% Inputs:
%  Coordinates = nodal coordinates [x y]
%  Contour = nodal solution field [c]
%  NodesOnElement = element connectivity [n1 n2 ... nnen]
%  SurfacesI = interface connectivity [n1R n2R n1L n2L eL eR edgeL edge R]
%              values may contain -1's for first 4 columns of SurfacesI,
%              then assumes the entire edge is the segment
%  numSI = number of interface segments
%  {PlotID} = [ModelID {clfyn subp numsubp}]
%      ModelID = figure ID for plot window (1)
%      clfyn = 1 for clearing figure before plot, retain otherwise (1)
%      subp = subplot ID for current contour plot (1)
%      numsubp = number of subplots for figure window (1)
%  {elemlist} = list of specific segments to plot [e1 e2 ... en]' (1:numSI)
%  {MNE_ID} = [meshyn nodeyn elemyn] (1 0 0)
%      meshyn = 1 for plotting mesh lines, no mesh lines otherwise
%      nodeyn = 1 for plotting node ID on figure, no ID on plot otherwise
%      elemyn = 1 for plotting elem ID on figure, no ID on plot otherwise
%  {zcont3D} = 1 for using Contour as the z-coordinate for the plot,
%              creates a flat contour plot otherwise (0)
%  {varargin} = {types},{props1},{props2},...,{propsn} 
%               pairs of commands for figure, axes, etc.
%      types = {'fig','node','mesh','cb','xlab','axes',...} is a list of
%              all types of property pairs that follow, in the correct
%              order
%      props1 = {{'PropertyName1',PropertyValue1,'PropertyName2',... is a
%              of all properties assigned to a given type
%      Examples: varargin = {'cb'},{'FontSize',16}
%                varargin = {'title'},{'String','Make a visible title'}
%                varargin = {'axes'},{'CLim',[-10 10]} % resets caxis
%                varargin = {'node'},{'FontSize',32} % sets node ID font
%                varargin = {'line'},{'LineWidth',2} % changes mesh lines
%                varargin = {'xlab'},{'FontSize',18} % changes xlabel
%                varargin = {'xlab','node'},{'FontSize',32},{'FontSize',18}
%                varargin = 'xlab','FontSize',32 % one type doesn't need {}

% Default arguments
switch nargin
    case {1,2,3,4}
        error('Must supply Coordinates, Contour, NodesOnElement, SurfacesI, and numSI')
    case 5
        ModelID = 1;
        clfyn = 1;
        subp = 1;
        numsubp = 1;
        numel = size(NodesOnElement,1);
        elemlist = (1:numSI)';
        meshyn = 1;
        nodeyn = 0;
        elemyn = 0;
        zcont3D = 0;
    case 6
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 1;
            subp = 1;
            numsubp = 1;
        end
        numel = size(NodesOnElement,1);
        elemlist = (1:numSI)';
        meshyn = 1;
        nodeyn = 0;
        elemyn = 0;
        zcont3D = 0;
    case 7
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 1;
            subp = 1;
            numsubp = 1;
        end
        numel = size(NodesOnElement,1);
        meshyn = 1;
        nodeyn = 0;
        elemyn = 0;
        zcont3D = 0;
    case 8
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 1;
            subp = 1;
            numsubp = 1;
        end
        numel = size(NodesOnElement,1);
        if length(MNE_ID) >= 3
            meshyn = MNE_ID(1);
            nodeyn = MNE_ID(2);
            elemyn = MNE_ID(3);
        else
            meshyn = 1;
            nodeyn = 0;
            elemyn = 0;
        end
        zcont3D = 0;
    otherwise
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 1;
            subp = 1;
            numsubp = 1;
        end
        numel = size(NodesOnElement,1);
        if length(MNE_ID) >= 3
            meshyn = MNE_ID(1);
            nodeyn = MNE_ID(2);
            elemyn = MNE_ID(3);
        else
            meshyn = 1;
            nodeyn = 0;
            elemyn = 0;
        end
end

%Select and clear figure
figure(ModelID);
% set(ModelID,'Position',[360   425   560   273])
% set(ModelID,'Position',[360   278   560   420]) %old HP
set(ModelID,'Position',[403   246   560   420])
if clfyn
clf(ModelID)
end

% Set nen = maximum number of nodes per element
nen = size(NodesOnElement,2)-1;

% Set nen = maximum number of nodes per element
if ~exist('MNE_ID','var')
    MNE_ID = [meshyn nodeyn elemyn];
end
        
% Open subplot
subplot(numsubp,1,subp)

colormap jet

% Set up arrays for node and element handles
if elemyn
    helem = zeros(length(elemlist),1);
else
    helem = 0;
end

% Plot mesh skeleton
plotMesh2(Coordinates,NodesOnElement,ModelID,1:numel-numSI,[MNE_ID 0 0 0],varargin)
axis on

% Nodes for edges
nloop3 = zeros(3,2);
nloop4 = zeros(4,2);
for i=1:3
    nloop3(i,1) = i;
    nloop3(i,2) = i+1;
end
nloop3(3,2) = 1;
for i=1:4
    nloop4(i,1) = i;
    nloop4(i,2) = i+1;
end
nloop4(4,2) = 1;

% Plot taus on interior edges
for elemp = 1:length(elemlist)
    
    ii_inter = elemlist(elemp);
    
    if sign(SurfacesI(ii_inter,1))>0 && sign(SurfacesI(ii_inter,2))>0
    xy = Coordinates(SurfacesI(ii_inter,1:2),:);   
    else
        elem = SurfacesI(ii_inter,5);
        edge = SurfacesI(ii_inter,7);
        nel = nnz(NodesOnElement(elem,1:nen));
        if (nel==3||nel==6)
        nodeA = NodesOnElement(elem,nloop3(edge,1));
        nodeB = NodesOnElement(elem,nloop3(edge,2));
        else
        nodeA = NodesOnElement(elem,nloop4(edge,1));
        nodeB = NodesOnElement(elem,nloop4(edge,2));
        end
    xy = Coordinates([nodeA,nodeB],:);   
    end
    c = ones(2,1)*Contour(ii_inter+numel-numSI);%ep_12(ii_inter);% % get value of tau for the edge
    if zcont3D == 1
    z = ones(2,1)*c;
    else
    z = zeros(2,1);
    end
    helem(ii_inter) = color_line3(xy(:,1)',xy(:,2)',z,c,'LineWidth',2);
    
end

% Final properties
shading interp
hx = xlabel('x');
hy = ylabel('y');
hz = zlabel('');
ht = title('');

axis equal
% axis([0 10 -3 3])
% axis off

hcb = colorbar('EastOutside','FontSize',13,'FontName','Arial','LineWidth',1);
% colorbar('EastOutside','FontSize',13,'FontName','Arial','LineWidth',1,'YTick',1e4*[-5 -2.5 0 2.5 5],'YTickLabel',[-5 -2.5 0 2.5 5])

% Set additional figure properties
if ~isempty(varargin)
if iscell(varargin{1,1}) % Multiple types of objects specified
    
    ObjTypes = varargin{1,1}; % header with list of object property types
    numObj = min(length(ObjTypes),length(varargin)-1);
    
    for obj = 1:numObj
        
        objtype = ObjTypes{1,obj}; % get title of object type
        switch objtype
            case {'fig','Fig','Figure','figure'}
                figstuff = varargin{1,1+obj}; % get list of figure properties
                set(gcf,figstuff{:}) % set figure properties
            case {'mesh','Mesh','line','Line','Lines','lin','lines'}
                linestuff = varargin{1,1+obj}; % get list of line properties
                set(findobj(gca,'Type','line'),linestuff{:}) % set line properties
            case {'axes'}
                axesstuff = varargin{1,1+obj}; % get list of axes properties
                set(gca,axesstuff{:}) % set axes properties
            case {'node','Node','nodes','Nodes','nodetext','NodeText'}
                if nodeyn
                nodestuff = varargin{1,1+obj}; % get list of node text properties
                set(hnode(1:hnode_count),nodestuff{:}) % set node text properties
                end
            case {'elem','Elem','elems','Elems','elements','Elements','elemtext','ElemText'}
                if elemyn
                elemstuff = varargin{1,1+obj}; % get list of elem text properties
                set(helem,elemstuff{:}) % set elem text properties
                end
            case {'xlabel','Xlabel','XLabel','xlab','Xlab'}
                xlabstuff = varargin{1,1+obj}; % get list of xlabel properties
                set(hx,xlabstuff{:}) % set xlabel properties
            case {'ylabel','Ylabel','YLabel','ylab','Ylab'}
                ylabstuff = varargin{1,1+obj}; % get list of ylabel properties
                set(hy,ylabstuff{:}) % set ylabel properties
            case {'zlabel','Zlabel','ZLabel','zlab','Zlab'}
                zlabstuff = varargin{1,1+obj}; % get list of zlabel properties
                set(hz,zlabstuff{:}) % set xlabel properties
            case {'title','Title'}
                titlestuff = varargin{1,1+obj}; % get list of title properties
                set(ht,titlestuff{:}) % set title properties
            case {'colorbar','CB','Colorbar','ColorBar','cb'}
                cbstuff = varargin{1,1+obj}; % get list of title properties
                set(hcb,cbstuff{:}) % set title properties
        end
        
    end
    
else % Single object only
        
        objtype = varargin{1,1}; % get title of object type
        switch objtype
            case {'fig','Fig','Figure','figure'}
                figstuff = varargin(1,2:end); % get list of figure properties
                set(gcf,figstuff{:}) % set figure properties
            case {'mesh','Mesh','line','Line','Lines','lin','lines'}
                linestuff = varargin(1,2:end); % get list of line properties
                set(findobj(gca,'Type','line'),linestuff{:}) % set line properties
            case {'axes'}
                axesstuff = varargin(1,2:end); % get list of axes properties
                set(gca,axesstuff{:}) % set axes properties
            case {'node','Node','nodes','Nodes','nodetext','NodeText'}
                if nodeyn
                nodestuff = varargin(1,2:end); % get list of node text properties
                set(hnode(1:hnode_count),nodestuff{:}) % set node text properties
                end
            case {'elem','Elem','elems','Elems','elements','Elements','elemtext','ElemText'}
                if elemyn
                elemstuff = varargin(1,2:end); % get list of elem text properties
                set(helem,elemstuff{:}) % set elem text properties
                end
            case {'xlabel','Xlabel','XLabel','xlab','Xlab'}
                xlabstuff = varargin(1,2:end); % get list of xlabel properties
                set(hx,xlabstuff{:}) % set xlabel properties
            case {'ylabel','Ylabel','YLabel','ylab','Ylab'}
                ylabstuff = varargin(1,2:end); % get list of ylabel properties
                set(hy,ylabstuff{:}) % set ylabel properties
            case {'zlabel','Zlabel','ZLabel','zlab','Zlab'}
                zlabstuff = varargin(1,2:end); % get list of zlabel properties
                set(hz,zlabstuff{:}) % set xlabel properties
            case {'title','Title'}
                titlestuff = varargin(1,2:end); % get list of title properties
                set(ht,titlestuff{:}) % set title properties
            case {'colorbar','CB','Colorbar','ColorBar','cb'}
                cbstuff = varargin(1,2:end); % get list of title properties
                set(hcb,cbstuff{:}) % set title properties
        end
    
end
end

% In case many numbers appear on colorbar legend:
set(ModelID, 'renderer', 'zbuffer');