% 04/30/2013

lEGDOFTa  = length(EGDOFTa );
lEGDOFTi  = length(EGDOFTi );
Ksdd(sinddd+1:sinddd+lEGDOFTa^2,:) = combvec(EGDOFTa,EGDOFTa)';
sinddd = sinddd + lEGDOFTa^2;
Ksdf(sinddf+1:sinddf+lEGDOFTa*lEGDOFTi,:) = combvec(EGDOFTa,EGDOFTi)';
sinddf = sinddf + lEGDOFTa*lEGDOFTi;
Ksff(sindff+1:sindff+lEGDOFTi^2,:) = combvec(EGDOFTi,EGDOFTi)';
sindff = sindff + lEGDOFTi^2;

% for locind1 = 1:nst
%     grow = EDOFT(locind1);
%     if grow <= neq && grow > 0
%         for locind2 = 1:nst
%             gcol = EDOFT(locind2);
%             if gcol <= neq && gcol > 0
%                 sinddd = sinddd + 1;
%                 Ksdd(sinddd,:) = [grow gcol]; %#ok<SAGROW>
%             elseif gcol > neq
%                 sinddf = sinddf + 1;
%                 Ksdf(sinddf,:) = [grow gcol-neq]; %#ok<SAGROW>
%             end %if gcol
%         end %ie2
%     elseif grow > neq
%         for locind2 = 1:nst
%             gcol = EDOFT(locind2);
%             if gcol > neq
%                 sindff = sindff + 1;
%                 Ksff(sindff,:) = [grow-neq gcol-neq]; %#ok<SAGROW>
%             end %if gcol
%         end %ie2
%     end %if grow
% end %ie1