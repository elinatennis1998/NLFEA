% Assemble Element Energy (Error) into Model Energy (Error)
%
% Tim Truster
% 7/2009
% UIUC

% for ener = 1:numEn
%     Energy(ener) = Energy(ener) + ElemE(ener); %#ok<AGROW>
% end
Energy(1:numEn) = Energy(1:numEn) + ElemE(1:numEn);