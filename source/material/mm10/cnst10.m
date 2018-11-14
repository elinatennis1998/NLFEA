% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine cnst10                            *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 3/22/12                     *
% c     *                                                              *
% c     *    gets the consistent tangent for a block                   *
% c     *                                                              *
% c     *                                                              *
% c     ****************************************************************
% c
function local_work...
    = cnst10(gp, span, history_np1, local_work)
%       use segmental_curves, only: max_seg_points
%       implicit integer (a-z)
% $add param_def
% $add include_tan_ek
%       integer, intent(in) :: iout
%       integer :: span, ncrystals(mxvl), hist_sz, gp
%       double precision :: history_n(span,hist_sz)
%       double precision :: history_np1(span,hist_sz)
% c
%       integer :: i
%
for i = 1:span;
    local_work.cep(i,:,:) = reshape(history_np1(i,1:36),6,6);
    local_work.cep(i,:,:) = local_work.cep(i,:,:)*...
        local_work.det_jac_block(i,gp)*local_work.weights(gp);
end
% c
%       return
end