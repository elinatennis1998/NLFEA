% c
% c
% c     ORNL:
% c
% c           Setup ORNL hardening
function [props, np1, n] = mm10_setup_ornl(props, np1, n)


% c     
% increment the total time
    time = n.u(1) + np1.tinc;
    np1.u(1) = time;
%       return
end