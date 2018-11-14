% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_set_e_nu                     *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 7/11/12                     *
% c     *                                                              *
% c     *    An annoying extra helper to enable the linear stiffness   *
% c     *    routines to play nicely with the CP model                 *
% c     *                                                              *
% c     *                                                              *
% c     ****************************************************************
% c
function [matnum,elnum,e,nu] = mm10_set_e_nu(matnum,elnum,e,nu)
%       use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
%       use crystal_data, only : c_array, angle_input, crystal_input,
%      &                              data_offset
% c
%       implicit integer (a-z)
% $add param_def
% c
%       integer :: matnum, elnum
%       real, intent(out) :: e, nu
% c
% c     Local
%       integer :: c, cnum, ncry, osn
% 
       e = 0.0;
       nu = 0.0;
       ncry = imatprp(101,matnum);
      for c = 1:ncry
% c                 Get the local crystal number
            if (imatprp(104,matnum) == 1), then
                  cnum = imatprp(105,matnum);
            elseif (imatprp(104,matnum) == 2), then
                  osn = data_offset(elnum); 
                  cnum = crystal_input(osn,c);
% c                 Couldn't do this earlier, so check here
%                   if ((cnum > max_crystals) | (cnum < 0)), then
%                    write (*,'("Crystal ", i3, " not valid")')cnum
%                         call die_gracefully
%                   end    % commented out since it is a write function
%             else
%                    write(*,*) "INSANITY IN mm10" 
%                    call die_gracefully
             end % partially commented out since it has a write function
                  e = e + (c_array(cnum).e);
                  nu = nu+ (c_array(cnum).nu);
      end
      e = e / ncry;
      nu = nu / ncry;
%       return
end