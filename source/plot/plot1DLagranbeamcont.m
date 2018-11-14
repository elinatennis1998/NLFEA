function helem = plot1DLagranbeamcont(P, Cont, nel, knots, numu, elem, helem)


SPL2D = zeros(numu,1);
CPL2D = zeros(numu,1);

%Generate list of curve points
ind = 0;
for i = 1:numu
    ind = ind + 1;
    r = knots(ind,1);
    shp = shl1d(r,1,nel-1);
    [shpw,shpt] = shp1dh(r,P(2,1)-P(1,1));
    SPL2D(i) = P'*shp';
    CPL2D(i) = Cont(:,1)'*shpw(4,:)' + Cont(:,2)'*shpt(4,:)';
end

%Plot the contour
hold on
helem(elem) = plot(SPL2D,CPL2D);
hold off