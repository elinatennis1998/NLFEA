function helem = plot2DLagransurfstre(P, Cont, nel, knots, numu, gridlin, ID, mateprop, compon,nelP,iel, zcont3D, elem, helem)
%(P, Cont, unum, vnum, numelU, numelV, pU, pV, numu, numv, Ulist, Vlist, gridlin, ID)
% plot3DNURBSsurf(Pw, U, V, nU, nV, pU, pV, numu, numv, facecol, gridlin)
%
% Tim Truster
% CEE Graduate Student
% UIUC
% 06/22/2009
%
%

SPL2D = zeros(numu,numu,2);
CPL2D = zeros(numu,numu);

%Generate list of curve points
ind = 0;
for j = 1:numu
    for i = 1:numu
        ind = ind + 1;
        r = knots(ind,1);
        s = knots(ind,2);
        [SurPoint,CPL2D(i,j)] = LagrSurfacePointStre(P,Cont,nel,r,s,2,mateprop,compon,nelP,iel);
        for dir = 1:2
            SPL2D(i,j,dir) = SurPoint(dir);
        end
    end
end

%Plot the list of curve points, control points, and knots
hold on
if zcont3D
surf(SPL2D(:,:,1),SPL2D(:,:,2),CPL2D,CPL2D,'EdgeColor',gridlin)%,zeros(numu,numu),CPL2D,'EdgeColor',gridlin)%
else
surf(SPL2D(:,:,1),SPL2D(:,:,2),zeros(numu,numu),CPL2D,'EdgeColor',gridlin)%,CPL2D,CPL2D,'EdgeColor',gridlin)%
end
uinc = round((numu-1)/2)+1;
vinc = round((numu-1)/2)+1;
if ID > 0
    if zcont3D
    helem(elem) = text(SPL2D(uinc,vinc,1),SPL2D(uinc,vinc,2),CPL2D(uinc,vinc),num2str(ID),'HorizontalAlignment','center');
    else
    helem(elem) = text(SPL2D(uinc,vinc,1),SPL2D(uinc,vinc,2),0,num2str(ID),'HorizontalAlignment','center');
    end
end
hold off