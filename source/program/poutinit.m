
    % Initialize lists for processing stresses/errors
    if SHist
        if ~exist('SHList','var')
            SHList = 1:stepmax;
        end
        if datastep > stepmax
        StreList = zeros(npstr,numnp,datastep);
        else
        StreList = zeros(npstr,numnp,length(SHList));
        end
        % start with first step that the user desires quantities
        stepS = 1;
    end
    if SEHist
        if ~exist('SEHList','var')
            SEHList = 1:stepmax;
        end
        if datastep > stepmax
        StreListE = zeros(nestr,numel,datastep);
        else
        StreListE = zeros(nestr,numel,length(SEHList));
        end
        % start with first step that the user desires quantities
        stepSE = 1;
    end
    if expliciterr
        if ~exist('EEHList','var')
            EEHList = stepmax;
        end
        EEList = zeros(numEn,length(EEHList));
        IeffList =zeros(5,numel,length(EEHList));
        % start with first step that the user desires quantities
        stepEE = 1;
    end
    if impliciterr
        if ~exist('IEHList','var')
            IEHList = stepmax;
        end
%         IEList = zeros(***,length(EEHList)); % implicit error is
%         specially implemented
        % start with first step that the user desires quantities
        stepIE = 1;
    end
    if PHist == 1
    %     StreList = zeros(2*(3*ndf-3),numel*nel,stepmax+1);
        PlasList = zeros(12,numel*nel,datastep);
    end
    if EHist == 1
        EnerList = zeros(1,datastep);
    end
    if PDHist == 1
        PDList = zeros(1,datastep);
    end
    if SSHist == 1
        SSList = zeros(13,datastep);
    end
    if IFHist == 1
        isw = 61;
        FormIData
    end
    if JHist == 1
        JintList = zeros(ndm,JDImax,datastep);
    end
    