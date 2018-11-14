% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_calc_grads                   *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                             2d version                       *
% c     *                   last modified: 08/14/14                    *
% c     *                                                              *
% c     *    calculate the gradient of Re.T (Fe) through a linear      *
% c     *    fit                                                       *
% c     *                                                              *
% c     ****************************************************************
% c
function gradFes =...
    mm10_calc_grads2(ngp, elem_type, order, geonl, rot_blk, jac,...
    Rps,nel)
%       implicit none
%       integer :: ngp, elem_type, order
%       double precision, dimension(9,ngp) :: rot_blk, Rps
%       double precision, dimension(27,ngp) :: gradFes
%       double precision, dimension(3,3) :: jac
%       logical :: geonl
% c
%       integer :: i, a, b, lwork, info
Rt = zeros(ngp,3,3); 
%       double precision, dimension(3,3) :: jacinv
intermat = zeros(ngp,3); 
%       double precision :: weight, fact
%       double precision, dimension(8) :: work
grads = zeros(3,3,3); 
% c
% c
% c           Get R components and stick in the right place
      if ngp == 1 % squeeze function outside makes the array dimensioned wrong
          rot_blk = rot_blk';
      end
      if (geonl)
        jacinv = jac;
        jacinv = inv(jacinv);
        for i = 1:ngp
          Rt(i,1:3,1:3) = reshape(Rps(1:9,i),3,3)*...
              reshape(rot_blk(1:9,i),3,3)';
        end
      else
        for i=1:ngp
          Rt(i,1:3,1:3) = reshape(Rps(1:9,i), 3,3);
        end
      end
% c
% c     For each Rt component create an interpolation, solve for the
% c     coefficients, and store the gradient
      for a=1:2
        for b=1:2
%            intermat = 0.0;
%            RHS = 0.0;
RHS = zeros(ngp,1);
          for i=1:ngp
% c           1-3 are the coordinates
              if nel == 3 || nel == 6
              [~,r,s] =  intpntt(i,ngp,0);
              else
              [~,r,s] =  intpntq(i,ngp,0);
              end
              intermat(i,1:2) = [r s];
            intermat(i,3) = 1.0;
            RHS(i) = Rt(i,a,b);
          end
% c           Solve with LAPACK
%           lwork = 8;
          RHS = intermat\RHS;
%           ['N',  ~, 4, 1, ~, ~, RHS, ngp, work, ~,...
%               info] = DGELS('N',  ngp, 4, 1, intermat, ngp, RHS,...
%               ngp, work, lwork, info); % http://tinyurl.com/q24e695 describes DGELS
% c           Extract coefs
%           if (info ~= 0)
%             write (*,*) "Error finding least squares mm10 grad."
%           end
% c           Get the gradient
          grads(a,b,1) = RHS(1);
          grads(a,b,2) = RHS(2);
% c           Take to the current coordinates
          if (geonl)
            grads(a,b,1:3) = jacinv*squeeze(grads(a,b,1:3));
          end
        end
      end
% c
% c     Flatten and store
      for i=1:ngp
        gradFes(1:27, i) = reshape(grads(1:3,1:3,1:3), 27, 1);
      end
% c
%       return
end