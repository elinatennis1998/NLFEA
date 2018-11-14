function shp = shltCSTT(ss)


shp = zeros(4,1);

shp(1)=ss(1);
shp(2)=ss(2);
shp(3)=1.d0-ss(1)-ss(2)-ss(3);
shp(4)=ss(3);
