classdef local_work_lin_ek < handle % makes it act like pass-by-reference
    % From include_lin_ek
    properties
       cp_stiff = zeros(1,6,6,1);
       cp_g_rot = zeros(1,3,3,1);
       cep = zeros(1,6,6);
       det_jac_block = zeros(1,1);
       weights = zeros(1,1);
       ncrystals = 1;
    end
   
end % classdef