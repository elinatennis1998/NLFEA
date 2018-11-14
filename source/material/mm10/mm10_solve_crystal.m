% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_solve_crystal                *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Advance a crystal from n to np1, store tangent, and      *
% c     *     store other output                                       *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, cut] = mm10_solve_crystal(props,...
    np1, n, cut, iout, fat,subcycle)
%       use mm10_defs
[~,max_uhard] = maxparamsCP;
vec1 = zeros(max_uhard,1);
vec2 = zeros(max_uhard,1);
asymmetric_assembly = 0;
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       logical :: cut, fat
%       integer :: iout
% c
      [props, np1, n, cut, J] = mm10_solve_strup(props, np1, n, cut,subcycle);
%  c      
      if (cut)
%         write(iout,*) "mm10 stress update failed"
        return
      end %if
% c
      [props, np1, n] = mm10_tangent(props, np1, n, J);
% c
% c      call mm10_ur_tangent(props, np1, n)
% c
% c      call mm10_num_tangent(props, np1, n)
% c
      if (~ asymmetric_assembly && ~ fat)
        np1.tangent = 0.5*(np1.tangent + transpose(np1.tangent));
      end
 [props, np1, n, np1.stress, np1.tau_tilde, vec1,vec2] = mm10_formvecs(props,np1,n,np1.stress,... 
 np1.tau_tilde,vec1,vec2);
% c
      [props, np1, n, vec1,vec2] = mm10_update_rotation(props, np1, n, vec1,vec2);
% c
      [props, np1, n, ~,~] = mm10_output(props, np1, n, vec1,vec2);
% c
%       return
end