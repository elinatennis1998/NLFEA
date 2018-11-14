% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_solve_strup                  *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Solve the stress update adaptively (if required)         *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, fail, J] = mm10_solve_strup(props, np1, n, fail,subcycle)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       logical :: fail
% c
%       type(crystal_state) :: curr
%       double precision, dimension(6) :: stress, ostress
%       double precision :: tt, ott
%       integer :: cuts, mcuts, i
%       double precision :: frac, step, mult
%       logical :: debug
 mcuts = 8; % parameter
 mult = 0.5;1/3;0.25; % parameter
 debug = 0; % debug = .false. % parameter
       cuts = 0;
% c
      oldversion = 1;0;
% c
      stress = n.stress;
%       tt = n.tau_tilde;
      ostress = stress;
%       ott = tt;
% c
      [props, np1, n] = mm10_setup(props, np1, n,subcycle);
      tt = n.tau_tilde(1:props.num_hard); %%% Mark modifies n.tau_tilde in setup for initial load step
      np1.tau_tilde = np1.tau_tilde(1:props.num_hard); % resize to actual number of hardening variables
      ott = tt;
      fail = 0; % fail = .false.
       
      if oldversion
% c
       frac = 0.0;
       step = 1.0;
% 
      while (frac < 1.0)
        curr =...
            mm10_setup_np1(reshape(np1.R, 9,1), np1.D*(step+frac),...
            np1.tinc*(step+frac), (np1.temp-n.temp)*(step+frac)...
            +n.temp, np1.step, np1.elem, np1.gp);
        [props, curr, n] = mm10_setup(props, curr, n,subcycle);
        curr.tau_tilde = curr.tau_tilde(1:props.num_hard); % resize to actual number of hardening variables
        [props, curr, n, stress, tt, fail, J] = ...
            mm10_solve(props, curr, n, stress, tt, fail, cuts, step);
        if (fail)
          if (debug)
%               write(*,*) "Adapting" % commented out write statement
          end % if
          stress = ostress;
          tt = ott;
          step = step * mult;
          cuts = cuts + 1;
            if (cuts > mcuts), break
            end % if
           fail = 0; % fail = .false.
        else
          ostress = stress;
          ott = tt;
          frac = frac + step;
        end %if
      end % while
% c
% c
      if (fail || any(isnan(tt)) || any(isnan(stress)))
%       write(props.out,*)" >>> Warning: mm10 implicit solution failed."
      fail = 1; % fail = .true.
      np1.stress = n.stress;
      np1.tau_tilde = n.tau_tilde;
%       return
      end
% c      
      np1.stress = stress;
      np1.tau_tilde = tt;
      np1.tt_rate = curr.tt_rate; % store rates of hardening also
      
      % Also, output the Jacobian J so that the tangent calculation can use
      % it without recomputing it.
      % Also allows numerical version computed by fsolve.
      
      else
          % do sub-incrementation like is used by kbc in mm05
          nsubinc = 1;
          notdone = 1;
          nol = n;
          np2 = np1;
          
          while notdone
              % Perform time integration for load step using nsubinc
              [props, np1, n, stress, tt, fail, J] = ...
                  mm10_solve_nsteps(nsubinc, props, np1, n, subcycle, stress, tt, fail, cuts);
              
              % Check if integrated stress/hardening is accurate/converged
            if (fail) || cuts <= -1
              stress = ostress;
              tt = ott;
              nsubinc = nsubinc / mult;
              cuts = cuts + 1;
                if (cuts > mcuts), break
                end % if
               fail = 0; % fail = .false.
              notdone = 1;
            else
              ostress = stress;
              ott = tt;
              notdone = 0;
            end %if
              
          end
          
          % Copy stuff
          
      end
      
% c
%       return
end