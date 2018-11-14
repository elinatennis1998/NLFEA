% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Setup voche law hardening      
function [props, np1, n] = mm10_setup_voche(props, np1, n)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c     
%       % No setup actually required, but define a mu_harden at state np1
%       % for the sake of the CP model
      np1.mu_harden = props.stiffness(6,6);
% c     
%       Return
end