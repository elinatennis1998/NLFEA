% Example batch file execution of NL_FEA_Program
% Tim Truster
% 10/06/2013
%
% Example problem is Cook's Membrane from the journal paper:
% http://dx.doi.org/10.1016/j.cma.2013.08.010
% 
% Description: The loading is applied in 4 equal increments. The input file
% with NCR=1 only applies two steps. The continuation file continues from
% the 2nd step to simulate the 3rd and 4th steps. The first restart file
% returns to step 2 and resimulates the last 2 steps. The second restart
% file returns to step 1, overrides the warning that history variables are
% needed, and simulates the last 3 steps over. Each time, the same value
% for y-displacement at node 2 (0.297394470862267) is calculated.

clear;
clc;

batchinter = 'batch';
NCR = 1;
batchname = 'Batch_Cooks_Mem_I.m';

NL_FEA_Program
DispList(2,2,2)

NCR = 2;
batchname = 'Batch_Cooks_Mem_C.m';

NL_FEA_Program
DispList(2,2,4)

NCR = 3;
batchname = 'Batch_Cooks_Mem_R1.m';

NL_FEA_Program
DispList(2,2,4)

NCR = 3;
batchname = 'Batch_Cooks_Mem_R2.m';

NL_FEA_Program
DispList(2,2,4)
