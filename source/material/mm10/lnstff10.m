% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine lnstff10                          *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 3/22/12                     *
% c     *                                                              *
% c     *   Get the linearized tangent for a block                     *
% c     *   TODO: include the elastic rotation                         *
% c     *                                                              *
% c     ****************************************************************
% c
function local_work = lnstff10(gp, span, local_work)
% $add param_def
% $add include_lin_ek
% I don't think local_work is available yet.
% c
%       integer :: gp, span
% c
%       integer :: i, ci, tc, a, b
%       double precision, dimension(6,6) :: totalC, Cci, Srot, Ct
%       double precision, dimension(3,3) :: g
% c
      for i = 1: span
       tc = 0;
       totalC = 0.0;
      for ci = 1: local_work.ncrystals(i)
            g = squeeze(local_work.cp_g_rot(i,1:3,1:3,ci));
            Ct = squeeze(local_work.cp_stiff(i,1:6,1:6,ci));
%             Srot = 0.0;
            Srot = mm10_RT2RVE(transpose(g));
            Cci = (Ct*transpose(Srot));
            Cci = (Srot*Cci);
            totalC = totalC + Cci;
            tc = tc + 1;
      end
      totalC = totalC / double(tc);
      for b = 1: 6
       for a = 1: 6
        local_work.cep(i,a,b) = totalC(a,b) * local_work.weights(gp) *...
            local_work.det_jac_block(i, gp);
       end
      end
      end
% 
% 
%       return
% 
end