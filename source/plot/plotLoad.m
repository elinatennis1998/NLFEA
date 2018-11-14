function plotLoad(Coordinates, FcMat, PlotID, maxarrow, varargin)
% Tim Truster
% 11/28/2013
% 
% Plot nodal loads computed from FormFELoad
% Optional arguments are in {}, defaults are in ().
%
% Default usage (for copy/paste):
% plotLoad(Coordinates,FcMat,1,0)
%
% Inputs:
%  Coordinates = nodal coordinates [x y]
%  FcMat = nodal loads; can be either a single column or multiple columns
%  {PlotID} = [ModelID {clfyn subp numsubp}]
%      ModelID = figure ID for plot window (1)
%      clfyn = 1 for clearing figure before plot, retain otherwise (0)
%      subp = subplot ID for current contour plot (1)
%      numsubp = number of subplots for figure window (1)
%               otherwise
%  {maxarrow} = physical length of arror of the maximum load in a 
%               coordinate direction (0=calculate)
%  {varargin} = {types},{props1},{props2},...,{propsn} 
%               pairs of commands for figure, axes, etc.
%      types = {'fig','node','mesh','cb','xlab','axes',...} is a list of
%              all types of property pairs that follow, in the correct
%              order
%      props1 = {{'PropertyName1',PropertyValue1,'PropertyName2',... is a
%              of all properties assigned to a given type
%      Examples: varargin = {'arr'},{'EdgeColor',[0 1 0]} % resets color of
%                           arrow edge lines
%                varargin = {'arr'},{'FaceColor',[0 1 0]} % resets color of
%                           arrow head
%                varargin = {'title'},{'String','Make a visible title'}
%                varargin = {'axes'},{'CLim',[-10 10]} % resets caxis
%                varargin = {'xlab'},{'FontSize',18} % changes xlabel
%                varargin = {'xlab','node'},{'FontSize',32},{'FontSize',18}
%                varargin = 'xlab','FontSize',32 % one type doesn't need {}


% Default arguments
switch nargin
    case 1
        error('Must supply Coordinates and FcMat')
    case 2
        ModelID = 1;
        clfyn = 0;
        subp = 1;
        numsubp = 1;
        maxarrow = 0;
    case 3
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 0;
            subp = 1;
            numsubp = 1;
        end
        maxarrow = 0;
    otherwise
        ModelID = PlotID(1);
        if length(PlotID) == 4
            clfyn = PlotID(2);
            subp = PlotID(3);
            numsubp = PlotID(4);
        else
            clfyn = 0;
            subp = 1;
            numsubp = 1;
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

% Open subplot
subplot(numsubp,1,subp)

ndf = size(FcMat,2);
ndm = min(size(Coordinates));
Coordinates = reshape(Coordinates,length(Coordinates),ndm);
numnp = size(Coordinates,1);
if ndf >= ndm
    ndf = ndm;
else
    ndf = 1;
end

% Set up arrays for handles
hload = zeros(length(Coordinates),ndf);
hload_count = zeros(1,ndf);

% Set up tolerance and arrow sizes
tol = 1e-12;
nodemax = max(Coordinates);
nodemin = min(Coordinates);
if maxarrow == 0
    loadmax = max(abs(FcMat));
    domainsize = norm(nodemax-nodemin);
    fraction = 0.1;
    maxarrow = fraction*domainsize;
    lennorm = maxarrow*ones(1,ndf)./loadmax(1:ndf);
else
    loadmax = max(abs(FcMat));
    lennorm = maxarrow*ones(1,ndf)./loadmax(1:ndf);
end
axis_limit = [nodemin-maxarrow*ones(1,ndm)
              nodemax+maxarrow*ones(1,ndm)];
axis_limit = reshape(axis_limit,1,2*ndm);
if ndm == 1
    axis2 = axis;
    axis_limit = [axis_limit axis2(3:4)];
end
axis(axis_limit)

headsize = 16;
I = eye(ndm);

x1 = zeros(1,max(2,ndm));
x2 = x1;



%%%% NOW: use quiver and quiver3

for node = 1:numnp
    
    if ndf == ndm
        loadnorm = abs(FcMat(node,:));
        for i = 1:ndm
            if loadnorm(i) > tol
                hload_count(i) = hload_count(i) + 1;
                x1(1:ndm) = Coordinates(node,1:ndm); % start at the actual node
                x2(1:ndm) = x1(1:ndm) + lennorm(i)*FcMat(node,i)*I(i,(1:ndm)); % arrow extends proportional to load
                headfact = lennorm(i)*loadnorm(i)/maxarrow; % scale headsize
                hload(hload_count(i),i) = arrow(x1,x2,headfact*headsize);
            end
        end
    else
        loadnorm = abs(FcMat(node,1));
        if loadnorm > tol
            hload_count = hload_count + 1;
            x1 = Coordinates(node,:); % start at the actual node
            x2 = x1 + lennorm*FcMat(node,1)*I(1,:); % arrow extends proportional to load
            headfact = lennorm*loadnorm/maxarrow; % scale headsize
            hload(hload_count) = arrow(x1,x2,headfact*headsize);
        end
    end

end

% Final properties
hx = xlabel('x');
hy = ylabel('y');
hz = zlabel('');
ht = title('');

axis equal
% axis([0 10 -3 3])
% axis off

% Reset color to red
for i = 1:ndf
set(hload(1:hload_count(i),i),'EdgeColor',[1 0 0],'FaceColor',[1 0 0]) % set dof marker properties
end

% Set additional figure properties
if ~isempty(varargin)
if iscell(varargin{1,1}) % Multiple types of objects specified
    
    ObjTypes = varargin{1,1}; % header with list of object property types
    numObj = min(length(ObjTypes),length(varargin)-1);
    
    for obj = 1:numObj
        
        objtype = ObjTypes{1,obj}; % get title of object type
        switch objtype
            case {'arrow','Arrow','Arr','arr'}
                arrstuff = varargin{1,1+obj}; % get list of dof marker properties
                for i = 1:ndf
                set(hload(1:hload_count(i),i),arrstuff{:}) % set dof marker properties
                end
            case {'fig','Fig','Figure','figure'}
                figstuff = varargin{1,1+obj}; % get list of figure properties
                set(gcf,figstuff{:}) % set figure properties
            case {'axes'}
                axesstuff = varargin{1,1+obj}; % get list of axes properties
                set(gca,axesstuff{:}) % set axes properties
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
        end
        
    end
    
else % Single object only
        
        objtype = varargin{1,1}; % get title of object type
        switch objtype
            case {'arrow','Arrow','Arr','arr'}
                arrstuff = varargin(1,2:end); % get list of dof marker properties
                for i = 1:ndf
                set(hload(1:hload_count(i),i),arrstuff{:}) % set dof marker properties
                end
            case {'fig','Fig','Figure','figure'}
                figstuff = varargin(1,2:end); % get list of figure properties
                set(gcf,figstuff{:}) % set figure properties
            case {'axes'}
                axesstuff = varargin(1,2:end); % get list of axes properties
                set(gca,axesstuff{:}) % set axes properties
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
        end
    
end
end

% In case many numbers appear on colorbar legend:
set(ModelID, 'renderer', 'zbuffer');
