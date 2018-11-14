% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formR1                       *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form R1                                                  *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, vec1, vec2, stress, tt, R1] = mm10_formR1(props, np1,...
     n, vec1, vec2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: R1
%       double precision :: tt
% c
%       double precision, dimension(6) :: dbarp
%       double precision, dimension(3) :: wp
%       double precision, dimension(6) :: symTW
% c
      [props, np1, n, vec1, vec2, stress, tt, dbarp] = mm10_form_dbarp(props,...
          np1, n, vec1, vec2, stress, tt);
      [props, np1, n, vec1, vec2, stress, tt, wp] = mm10_form_wp(props, np1,...
          n, vec1, vec2, stress, tt);
      work_vec1 = np1.D' - dbarp;
%
%	add R of creep strain due to pressure precipitation (Ran)
      if ( props.h_type==14 && abs(props.cp_031-1)<1.0e-5 )
          work_vec1 = mm10_halite_formRpp( props, work_vec1, stress, np1.tinc );
      end
%
      symTW = mm10_symSW(stress, wp);
% c
      R1 = transpose(stress) - n.stress' - props.stiffness*work_vec1...
          + 2.0 * symTW;
% c
%       return
end
