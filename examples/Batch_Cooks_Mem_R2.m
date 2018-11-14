% Example Restart file - intermediate step
% Tim Truster
% 10/06/2013
%
% Part of Batch_Cooks_Mem

NCR = 3;

currstep = 'n'; % y/n flag that tells the program to restart from last step in the restart data
rname = 'Restart0002.mat'; % filename of MATLAB .mat file, including path
stepin = 1; % step from which to restart the simulation
stepmax = 4; % maximum step number to continue simulation
hisover = 1; %  T/F flag which allows restart from an intermediate time level when hrmat DNE
mults = [0.25 0.5 0.75 1]'; % subsequent values of load factor lamda

% Other changes to default algorithmic parameters and material models and
% loads could also be applied in this file.