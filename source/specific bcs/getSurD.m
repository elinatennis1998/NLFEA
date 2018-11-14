function [SurfacesD] = getSurD(SurfacesL,Coordinates,D)
k=1;
for i = 1:size(SurfacesL,1)
   if Coordinates(SurfacesL(i,1),2)==D&&Coordinates(SurfacesL(i,2),2)==D
       SurfacesD(k,:) = SurfacesL(i,:);
       k= k+1;
   end
end