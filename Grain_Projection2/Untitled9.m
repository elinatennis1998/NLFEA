batchinter = 'batch';
batchname = 'Lel08_2d_Q4DG_Test';
% Allstuff = zeros(10*4+2,8);

bs = [1 2 3 4 6 12 24 48 96 96*2 96*4];
for jjj = 1:6%9%10:10%

bCrys = bs(jjj);
fprintf('starting file %s\n',bCrys)
FEA_Program
HillMandel2
Allstuff(:,jjj) = [interior boundary BoundEps(2,inters) BoundEps(3,inters) BoundEps(4,inters) BoundEps(5,inters) BoundEps(6,inters) BoundEps(7,inters) ...
     BoundEps(8,inters) BoundEps(9,inters) BoundEps(10,inters) BoundEps(11,inters)]';

end