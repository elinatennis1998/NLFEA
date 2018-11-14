function lint = IntPoint(nel)
% !
% !**********************************************************************
% 	Implicit None
% 	Integer Nel,Lint
        if nel == 3
            lint = 3;7;25; %minimum of 7 for all integrals in deformed state
        elseif nel == 4
            lint = 4;9;
%             lint = 100;16;
        elseif nel == 6
            lint = 13; %minimum of 13 for all integrals in deformed state
        else
%             lint = 9;
            lint = 25;
        end