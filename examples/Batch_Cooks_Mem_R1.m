% Example Restart file - restart from final step
% Tim Truster
% 10/06/2013
%
% Part of Batch_Cooks_Membrane

NCR = 3;

currstep = 'y'; % y/n flag that tells the program to restart from last step in the restart data
rname = 'Restart0002.mat'; % filename of MATLAB .mat file, including path
stepmax = 4; % maximum step number to continue simulation
mults = [0.25 0.5 0.75 1]'; % subsequent values of load factor lamda