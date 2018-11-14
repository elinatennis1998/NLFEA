function plot2DLagransurfebub(P, Cont, nel, knots, numu, gridlin, ID, zcont3D)

SPL2D = zeros(numu,numu,2);
CPL2D = zeros(numu,numu);

%Generate list of curve points
ind = 0;
for j = 1:numu
    for i = 1:numu
        ind = ind + 1;
        r = knots(ind,1);
        s = knots(ind,2);
        [SurPoint,CPL2D(i,j)] = LagrSurfacePointEbub(P,Cont,nel,r,s,2);
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
text(SPL2D(uinc,vinc,1),SPL2D(uinc,vinc,2),0,num2str(ID))
end
hold off