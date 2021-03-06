% c
% c
% c     MRR:
% c
% c           Form intermediate arrays for faster calculations
function [props, np1, n, stress, tt, arr1, arr2] = mm10_a_mrr(props, np1,...
    n, stress, tt, arr1, arr2,both)

        [props, np1, n, stress, tt, arr1(1:props.num_hard,1)] = mm10_dgdt_mrr(props, np1,...
            n, stress, tt);
        if both == 2
          [props, np1, n, stress, tt, arr2(1:props.nslip,1:props.num_hard)] = mm10_dgdh_mrr(props, np1,...
            n, stress, tt);
        end

end