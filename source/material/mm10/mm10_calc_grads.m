% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_calc_grads                   *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 06/22/15 (ON)               *
% c     *                                                              *
% c     *    calculate the gradient of Re.T (Fe) through a linear      *
% c     *    fit                                                       *
% c     *                                                              *
% c     ****************************************************************
% c
function gradFes =...
    mm10_calc_grads(ngp, elem_type, order, geonl, rot_blk, jac,...
    Rps,nel)

Rt = zeros(ngp,3,3);
intermat = zeros(ngp,4);
grads = zeros(3,3,3);

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
for a=1:3
    for b=1:3
        
        %         Loop:
        for i=1:ngp
            %             % c           1-3 are the coordinates
            %             [~,ss] =  intpntb(i,ngp,0);
            if nel == 4 || nel == 10
                [~,intermat(i,1:3)] =  int3d_t(i,ngp,0);
            else
                [~,intermat(i,1:3)] =  intpntb(i,ngp,0);
            end
        end
        
        %       Vector:
        intermat(:,4) = 1.0;
        LHS(:,1) = Rt(:,a,b);
        
        RHS = intermat\LHS;
        
        grads(a,b,1) = RHS(1);
        grads(a,b,2) = RHS(2);
        grads(a,b,3) = RHS(3);
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