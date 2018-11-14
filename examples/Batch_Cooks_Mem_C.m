% Example Continuation file
% Tim Truster
% 10/06/2013
%
% Part of Batch_Cooks_Mem

NCR = 2;

stepmax = 4; % maximum step number to continue simulation
mults = [mults' 0.75 1]; % subsequent values of load factor lamda

% Other changes to default algorithmic parameters and material models and
% loads could also be applied in this file.