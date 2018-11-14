% Tim Truster
% 01/25/2015

LineSerType = 2;
% lsOn = 1;

if Rnorm > Residtol*1e5 && lsOn

    if LineSerType == 1

    else

    % Line search using N-R

    Gratio = 1/2*10^-2; %epsilon in G(s) equation
    sLS = 0.9;0;1;
    jter = 0;
    denomG = (del_ModelDx'*Fd1); % (del_ModelDx'*Kn*del_ModelDx);
    ModelDx0 = ModelDx;
    
    % Get G
    ModelDx = ModelDx0 + sLS*del_ModelDx;
    isw = 6;
    FormFE
    Gval = del_ModelDx'*Fd1;
    
    Gtol = Gratio*abs(Gval);
    
    % iteration 1
    jter = jter + 1;
    del_s = Gval/denomG;
    sLS = sLS + del_s;
        
    if sLS <= 1
        
        % Get G
        ModelDx = ModelDx0 + sLS*del_ModelDx;
        isw = 6;
        FormFE
        Gval = del_ModelDx'*Fd1;
        
        while abs(Gval) > Gtol && jter < 5 && sLS <= 1
            jter = jter + 1;
            del_s = Gval/denomG;
            sLS = sLS + del_s;
            % Get G
            ModelDx = ModelDx0 + sLS*del_ModelDx;
            isw = 6;
            FormFE
            Gval = del_ModelDx'*Fd1;
        end
    
    end
    
    % Catch if sLS goes outside the allowable [0,1] range
    if sLS > 1 || sLS <= 0
        sLS = 1;
    else
        sLS;
    end
    
    ModelDx = ModelDx0;

    end
    
else
    
    sLS = 1;

end % Rnorm
