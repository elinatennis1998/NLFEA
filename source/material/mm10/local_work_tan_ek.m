classdef local_work_tan_ek < handle % makes it act like pass-by-reference
    % From include_lin_ek
    properties
       cp_stiff = zeros(1,6,6,8);
       cp_g_rot = zeros(1,3,3,8);
       cep = zeros(1,6,6);
       det_jac_block = zeros(1,8);
       weights = zeros(1,8);
       ncrystals = 1;
    end
   
end % classdef