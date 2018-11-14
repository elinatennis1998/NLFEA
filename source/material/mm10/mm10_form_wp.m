% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Form w_p (in the current configuration)
function [props, np1, n, vec1, vec2, stress, tt, w] = ...
   mm10_form_wp(props, np1, n, vec1, vec2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(3) :: w
%       double precision :: tt
% c
%       integer :: i
% c
%       double precision :: mm10_slipinc
% c
      if (props.h_type == 1) % voche
       wmat = np1.dg * (ones(3,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.qc;
       w = sum(wmat(:,1:props.nslip),2);
      elseif (props.h_type == 2) % MTS
      rss = stress*np1.ms;
      rss2 = rss.*sign(real(rss));
       wmat = np1.dg * (ones(3,1)*((rss2/tt).^(props.rate_n).*sign(real(rss)))).*np1.qc;
       w = sum(wmat(:,1:props.nslip),2);
      elseif (props.h_type == 3) % User
       wmat = np1.dg * (ones(3,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.qc;
       w = sum(wmat(:,1:props.nslip),2);
      elseif (props.h_type == 4) % ORNL
%         w2 = np1.qc*vec1(1:props.nslip);
        w = zeros(3,1);
        gam = vec1*0;
        for i = 1:props.nslip
            gam(i) = mm10_slipinc_ornl(props, np1, n, stress, tt, i);
            w = w + gam(i)*...
                  np1.qc(1:3,i);
        end
%         if norm(w2-w)>1e-16
%             w2;
%         else
%             w = w2;
%         end
      elseif (props.h_type == 5) % AADP
        w = zeros(3,1);
        for i = 1:props.nslip
            w = w + mm10_slipinc_aadp(props, np1, n, stress, tt, i)*...
                  np1.qc(1:3,i);
        end
elseif (props.h_type == 6) % MTS_Omar
    rss = stress*np1.ms-np1.backstress_omar;
    rss2 = rss.*sign(real(rss));
       wmat = np1.dg * (ones(3,1)*((rss2/tt).^(props.rate_n).*sign(real(rss)))).*np1.qc;
       w = sum(wmat(:,1:props.nslip),2);
      elseif (props.h_type == 7) % MRR
%         w2 = np1.qc*vec1(1:props.nslip);
        w = zeros(3,1);
        gam = vec1*0;
        for i = 1:props.nslip
            gam(i) = mm10_slipinc_mrr(props, np1, n, stress, tt, i);
            w = w + gam(i)*...
                  np1.qc(1:3,i);
        end
%         if norm(w2-w)>1e-16
%             w2;
%         else
%             w = w2;
%         end
      elseif (props.h_type == 8) % DYS
        w = zeros(3,1);
        for i = 1:props.nslip
            w = w + mm10_slipinc_dys(props, np1, n, stress, tt, i)*...
                  np1.qc(1:3,i);
        end
      elseif (props.h_type == 12) % Plugin Library
          w = mm10_plugin_lib(8,props,{props, np1, n, vec1, vec2, stress, tt});
      elseif (props.h_type == 9) % Ran_HCP
%         w2 = np1.qc*vec1(1:props.nslip);
        w = zeros(3,1);
        gam = vec1*0;
        for i = 1:props.nslip
            gam(i) = mm10_slipinc_DJGM(props, np1, n, stress, tt, i);
            w = w + gam(i)*...
                  np1.qc(1:3,i);
        end
      elseif (props.h_type == 14) % halite
%         w2 = np1.qc*vec1(1:props.nslip);
        w = zeros(3,1);
        gam = vec1*0;
        for i = 1:props.nslip
            gam(i) = mm10_slipinc_halite(props, np1, n, stress, tt, i);
            w = w + gam(i)*...
                  np1.qc(1:3,i);
        end
%         if norm(w2-w)>1e-16
%             w2;
%         else
%             w = w2;
%         end
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
%       return
end