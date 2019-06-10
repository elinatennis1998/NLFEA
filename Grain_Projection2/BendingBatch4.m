batchinter = 'batch';
% batchname = 'BendingTest9';
batchname = 'BendingTest4';
% Allstuff = zeros(10*4+2,8);

bs = 24;4;6;12;16;2;%[1 2 3 4 6 12 24 48 96 96*2 96*4];
for jjj = 1:1

bCrys = bs(jjj);
fprintf('starting file %s\n',bCrys)
xxs = zeros(16,3*bCrys);
xys = zeros(16,3*bCrys);
uxs = zeros(16,3*bCrys);
uys = zeros(16,3*bCrys);
sxxs = zeros(16,3*bCrys);
syys = zeros(16,3*bCrys);
sxys = zeros(16,3*bCrys);
exxs = zeros(16,3*bCrys);
eyys = zeros(16,3*bCrys);
exys = zeros(16,3*bCrys);
plotp = 0;
summm = 0;
sssss = 1;
FEA_Program
RegList = 1:nummatCG;
HillMandel4
% Allstuff(:,jjj) = [interior boundary BoundEps(2,inters) BoundEps(3,inters) BoundEps(4,inters) BoundEps(5,inters) BoundEps(6,inters) BoundEps(7,inters) ...
%      BoundEps(8,inters) BoundEps(9,inters) BoundEps(10,inters) BoundEps(11,inters)]';

end
% close all
% plotNodeCont2(Coordinates,StreList(1,:)',NodesOnElement,1,(1:size(NodesOnElement,1)),[1 0 0],0,[3 4 6 9 0],{'title'},{'String','Sxx'})
% plotNodeCont2(Coordinates,StreList(2,:)',NodesOnElement,2,(1:size(NodesOnElement,1)),[1 0 0],0,[3 4 6 9 0],{'title'},{'String','Syy'})
% plotNodeCont2(Coordinates,StreList(3,:)',NodesOnElement,3,(1:size(NodesOnElement,1)),[1 0 0],0,[3 4 6 9 0],{'title'},{'String','Sxy'})
% plotNodeCont2(Coordinates,Node_U_V(:,1),NodesOnElement,4,(1:size(NodesOnElement,1)),[1 0 0],0,[3 4 6 9 0],{'title'},{'String','Ux'})
% plotNodeCont2(Coordinates,Node_U_V(:,2),NodesOnElement,5,(1:size(NodesOnElement,1)),[1 0 0],0,[3 4 6 9 0],{'title'},{'String','Uy'})