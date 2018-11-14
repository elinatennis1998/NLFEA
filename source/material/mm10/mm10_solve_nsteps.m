% c
% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_solve_nsteps                 *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 08/23/15                    *
% c     *                                                              *
% c     *     Perform sub-incremental update of stresses               *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, fail, J] = mm10_solve_nsteps(nsubinc, props,...
    np1, n, subcycle, stress, tt, fail, cuts)
%
[~,max_uhard] = maxparamsCP;
vec1 = zeros(max_uhard,1);
vec2 = zeros(max_uhard,1);
% c
      ostress = stress;
      ott = tt;
      ott_rate = 0*np1.tt_rate;

      incsize = 1/nsubinc;
      prev = n; % copy values out of n to preserve it

      for sub = 1:nsubinc
        curr =...
            mm10_setup_np1(reshape(np1.R, 9,1), np1.D*incsize,...
            np1.tinc*incsize, (np1.temp-n.temp)*incsize...
            +prev.temp, np1.step, np1.elem, np1.gp);
        [props, curr, prev] = mm10_setup(props, curr, prev, subcycle);
        curr.tau_tilde = curr.tau_tilde(1:props.num_hard); % resize to actual number of hardening variables
        [props, curr, prev, stress, tt, fail, J] = ...
            mm10_solve(props, curr, prev, stress, tt, fail, cuts, incsize);
        if (fail)
          break
        else
          ostress = stress;
          ott = tt;
          curr.stress = stress; %transpose so that the vectorization works
          curr.tau_tilde = tt'; %transpose so that the vectorization works
          curr.Rp = prev.Rp;
%           [props, curr, prev, curr.stress, curr.tau_tilde, vec1,vec2] = mm10_formvecs(props,curr,prev,curr.stress, curr.tau_tilde, ... 
%             vec1,vec2);
% % c
%           [props, curr, prev, vec1,vec2] = mm10_update_rotation(props, curr, prev, vec1,vec2);
          prev = curr;
          ott_rate = ott_rate + incsize*curr.tt_rate; % store rates of hardening also
        end %if
      end % do
% c
% c
      if (fail || any(isnan(tt)) || any(isnan(stress)))
%       write(props.out,*)" >>> Warning: mm10 implicit solution failed."
      fail = 1; % fail = .true.
      np1.stress = n.stress;
      np1.tau_tilde = n.tau_tilde;
%       return
      else
% c      
      np1.stress = stress;
      np1.tau_tilde = tt;
      np1.tt_rate = ott_rate; % store rates of hardening also
      end
% c
%       return
end