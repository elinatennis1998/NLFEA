% c
% c
% c     MRR:
% c
% c           Setup MRR hardening
function [props, np1, n] = mm10_setup_mrr(props, np1, n)


% increment the total time
    time = n.u(1) + np1.tinc;
    np1.u(1) = time;
% c     
%       return
end