% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c      from mm10_form.f
% c      
% c           Calculate the resolved shear along system i     
function [mm10_rs] = ...
    mm10_rs(np1, stress, i)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       integer :: i
% c
%       double precision :: mm10_rs
% c
      mm10_rs = stress*np1.ms(1:6,i);
% c
%       return
end