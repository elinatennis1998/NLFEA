% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Form w_bar_p
function [props, np1, n, vec1, vec2, stress, tt, wbar] = ...
    mm10_form_wbarp(props, np1, n, vec1, vec2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(3) :: wbar
%       double precision :: tt
% c
%       integer :: i
% c
%       double precision :: mm10_slipinc
% c
      if (props.h_type == 1) % voche
       wbarmat = np1.dg * (ones(3,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.qs;
       wbar = sum(wbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 2) % MTS
       wbarmat = np1.dg * (ones(3,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.qs;
       wbar = sum(wbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 3) % User
       wbarmat = np1.dg * (ones(3,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.qs;
       wbar = sum(wbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 4) % ORNL
        wbar = np1.qs*vec1(1:props.nslip);
%         wbar = zeros(3,1);
%         for i = 1:props.nslip
%             wbar = wbar + mm10_slipinc_ornl(props, np1, n, stress, tt, i)*...
%                   np1.qs(1:3,i);
%         end
%         if norm(wbar2-wbar)>1e-9
%             wbar2;
%         else
%             wbar = wbar2;
%         end
      elseif (props.h_type == 5) % AADP
        wbar = 0;
        for i = 1:props.nslip
            wbar = wbar + mm10_slipinc_aadp(props, np1, n, stress, tt, i)*...
                  np1.qs(1:3,i);
        end
elseif (props.h_type == 6) % MTS_omar copy
    teff_omar = stress*np1.ms-np1.backstress_omar;
    wbarmat = np1.dg * (ones(3,1)*((teff_omar/tt).^(props.rate_n).*sign(teff_omar))).*np1.qs;
    wbar = sum(wbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 7) % MRR
        wbar = np1.qs*vec1(1:props.nslip);
%         wbar = zeros(3,1);
%         for i = 1:props.nslip
%             wbar = wbar + mm10_slipinc_mrr(props, np1, n, stress, tt, i)*...
%                   np1.qs(1:3,i);
%         end
%         if norm(wbar2-wbar)>1e-9
%             wbar2;
%         else
%             wbar = wbar2;
%         end
      elseif (props.h_type == 8) % DYS
        wbar = 0;
        for i = 1:props.nslip
            wbar = wbar + mm10_slipinc_dys(props, np1, n, stress, tt, i)*...
                  np1.qs(1:3,i);
        end
      elseif (props.h_type == 9) % Ran_HCP
        wbar = np1.qs*vec1(1:props.nslip);
      elseif (props.h_type == 12) % Plugin Library
          wbar = mm10_plugin_lib(7,props,{props, np1, n, vec1, vec2, stress, tt});
      elseif (props.h_type == 14) % halite
        wbar = np1.qs*vec1(1:props.nslip);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
%       return
end